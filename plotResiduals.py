#!/usr/bin/env python3
"""
Static plot of residuals from a completed simpleFoam log file.

Usage:
    python static_residual_plot.py --log path/to/log.simpleFoam --stride 100
"""

import argparse
import re
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt

# ─── Argument Parser ──────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="Static residual plot from simpleFoam log.")
    parser.add_argument(
        "--log", type=Path, default=Path("log.simpleFoam"),
        help="Path to log.simpleFoam file (default: log.simpleFoam)"
    )
    parser.add_argument(
        "--stride", type=int, default=1,
        help="Plot every N-th point to reduce load (default: 1 = no skipping)"
    )
    return parser.parse_args()

# ─── Residual Parser ──────────────────────────────────────────────────────────

def parse_residuals(lines):
    """
    Extract time and residuals from simpleFoam log lines.

    Returns:
        dict[str, tuple[list[float], list[float]]]: {field: ([time], [residual])}
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

# ─── Plotting Function ────────────────────────────────────────────────────────

def plot_static_residuals(data: dict, stride: int):
    plt.figure(figsize=(10, 6))
    plt.yscale("log")
    plt.xlabel("Time")
    plt.ylabel("Initial Residual")
    plt.title(f"Residuals (stride = {stride})")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)

    for idx, (field, (times, residuals)) in enumerate(sorted(data.items())):
        if len(times) == 0:
            continue
        plt.plot(
            times[::stride],
            residuals[::stride],
            label=field
        )

    plt.legend()
    plt.tight_layout()
    plt.show()

# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if not args.log.exists():
        print(f"[ERROR] Log file not found: {args.log}")
        return

    try:
        with args.log.open("r") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"[ERROR] Failed to read log file: {e}")
        return

    print(f"[INFO] Parsing residuals from: {args.log}")
    data = parse_residuals(lines)

    if not data:
        print("[WARNING] No residual data found.")
        return

    print(f"[INFO] Found fields: {', '.join(data.keys())}")
    plot_static_residuals(data, args.stride)

if __name__ == "__main__":
    main()
