#!/usr/bin/env python3
"""
Live plot residuals from a running simpleFoam simulation with optional stride control.

Usage:
    python live_residual_plot.py --log path/to/log.simpleFoam --stride 100
"""

import argparse
import time
import re
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ─────────────────────────────────────────────────────────────────────
# Argument Parsing
# ─────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="Live plot residuals from simpleFoam log.")
    parser.add_argument(
        "--log", type=Path, default=Path("log.simpleFoam"),
        help="Path to simpleFoam log file (default: log.simpleFoam)"
    )
    parser.add_argument(
        "--stride", type=int, default=1,
        help="Plot every N-th point (default: 1 = no skipping)"
    )
    return parser.parse_args()

# ─────────────────────────────────────────────────────────────────────
# Residual Parsing
# ─────────────────────────────────────────────────────────────────────

def parse_residuals(lines):
    """
    Extract time and residuals from lines of a simpleFoam log file.

    Returns:
        dict[str, tuple[list[float], list[float]]] - field name mapped to (time list, residual list)
    """
    data = defaultdict(lambda: ([], []))
    current_time = None

    time_re = re.compile(r"^Time = (\d+\.?\d*)")
    residual_re = re.compile(r"Solving for (\w+), Initial residual = ([\deE\+\-\.]+)")

    for line in lines:
        time_match = time_re.match(line)
        if time_match:
            current_time = float(time_match.group(1))
            continue

        residual_match = residual_re.search(line)
        if residual_match and current_time is not None:
            field, residual = residual_match.groups()
            data[field][0].append(current_time)
            data[field][1].append(float(residual))

    return data

# ─────────────────────────────────────────────────────────────────────
# Live Plotting Class
# ─────────────────────────────────────────────────────────────────────

class LiveResidualPlotter:
    def __init__(self, log_path: Path, stride: int):
        self.log_path = log_path
        self.offset = 0
        self.stride = max(1, stride)
        self.data = defaultdict(lambda: ([], []))
        self.colors = ["blue", "green", "red", "orange", "purple", "cyan", "magenta"]
        self.fig, self.ax = plt.subplots()

    def read_new_lines(self):
        try:
            with self.log_path.open("r") as f:
                f.seek(self.offset)
                lines = f.readlines()
                self.offset = f.tell()
                return lines
        except FileNotFoundError:
            print(f"[ERROR] Log file not found: {self.log_path}")
            return []

    def update(self, frame):
        new_lines = self.read_new_lines()
        if not new_lines:
            return

        new_data = parse_residuals(new_lines)

        for field, (times, residuals) in new_data.items():
            self.data[field][0].extend(times)
            self.data[field][1].extend(residuals)

        self.ax.clear()
        self.ax.set_yscale("log")
        self.ax.set_title(f"Live Residuals (stride={self.stride})")
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel("Initial Residual")
        self.ax.grid(True, which="both", linestyle="--", alpha=0.5)

        for idx, (field, (times, residuals)) in enumerate(self.data.items()):
            if not times:
                continue
            color = self.colors[idx % len(self.colors)]
            # Only plot every N-th point
            times_strided = times[::self.stride]
            residuals_strided = residuals[::self.stride]
            self.ax.plot(times_strided, residuals_strided, label=field, color=color)

        self.ax.legend(loc="upper right")

    def run(self):
        print(f"[INFO] Watching log file: {self.log_path}")
        print(f"[INFO] Plotting every {self.stride} points.")
        print(f"[INFO] Press Ctrl+C to stop.")
        ani = FuncAnimation(self.fig, self.update, interval=1000)
        plt.tight_layout()
        plt.show()

# ─────────────────────────────────────────────────────────────────────
# Main Script Execution
# ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    if not args.log.exists():
        print(f"[ERROR] Cannot find log file: {args.log}")
        return

    plotter = LiveResidualPlotter(args.log, stride=args.stride)
    try:
        plotter.run()
    except KeyboardInterrupt:
        print("\n[INFO] Interrupted by user — exiting.")
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")

if __name__ == "__main__":
    main()
