#!/usr/bin/env python3
"""
Post-process OpenFOAM outputs:
- Extract head losses, mass flow rates, and max yPlus from each simulation case.
- Assemble results into CSV tables.
- Optionally plot iteration data.

Directory structure:
  ./outputs/<caseName>/postProcessing
Outputs written to:
  ./postProcessingOutputs/
"""

import sys
import argparse
import re
from pathlib import Path
from typing import Optional, Tuple, List, Dict

# ────────────────────────────────
# Configuration Constants
# ────────────────────────────────
OUTPUTS_DIR = Path("outputs")
POSTPROC_DIR = Path("postProcessingOutputs")

# ────────────────────────────────
# Argument Parsing
# ────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract OpenFOAM post-processing data.")
    parser.add_argument("--write", action="store_true", help="Write results to CSV")
    parser.add_argument("--plot", action="store_true", help="Plot results")
    return parser.parse_args()

# ────────────────────────────────
# File Utilities
# ────────────────────────────────

def read_last_line(path: Path) -> str:
    """
    Read the last non-empty, non-comment line of a file.
    Optimized to avoid reading entire file unnecessarily.
    """
    with path.open("r") as f:
        for line in reversed(f.readlines()):
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                return stripped
    raise ValueError(f"No valid data in {path}")

def read_value_from_last_line(path: Path, index: int = -1) -> float:
    """
    Extract a specific token (by index) from the last valid line of a file.
    """
    line = read_last_line(path)
    tokens = line.split()
    if len(tokens) < abs(index) + 1:
        raise ValueError(f"Not enough tokens in {path}: {tokens}")
    return float(tokens[index])

def read_values(file_path: Path, iter_col: int = 0, value_col: int = -1) -> Tuple[List[float], List[float]]:
    """
    Read columns of iteration and value data from a file.
    Returns (iterations, values).
    """
    if not file_path.exists():
        raise FileNotFoundError(f"{file_path} not found")

    with file_path.open() as f:
        lines = f.readlines()

    valid_lines = [line.strip() for line in lines if line.strip() and not line.startswith("#")]
    data = [line.split() for line in valid_lines]

    iterations, values = [], []
    for line in data:
        try:
            iterations.append(float(line[iter_col]))
            values.append(float(line[value_col]))
        except (IndexError, ValueError):
            raise ValueError(f"Malformed line in {file_path}: {line}")

    return iterations, values

def safe_glob(base: Path, pattern: str) -> Path:
    """
    Find the first file that matches a pattern, or raise.
    """
    try:
        return next(base.glob(pattern))
    except StopIteration:
        raise FileNotFoundError(f"No file matching {pattern} in {base}")

def get_required_files(case_dir: Path) -> Dict[str, Path]:
    """
    Collect paths to required post-processing output files.
    """
    base = case_dir / "postProcessing"
    return {
        'solver': safe_glob(base, "solverInfo1/0/solverInfo.dat"),
        'p_in': safe_glob(base, "averageTotalPressureInlet1/0/surfaceFieldValue.dat"),
        'p_out': safe_glob(base, "averageTotalPressureOutlet1/0/surfaceFieldValue.dat"),
        'm_flow': safe_glob(base, "massFlowRateOutlet1/0/surfaceFieldValue.dat"),
        'y_plus': safe_glob(base, "yPlus1/0/yPlus.dat"),
    }

# ────────────────────────────────
# Data Extraction
# ────────────────────────────────

def extract_case_data(case_dir: Path) -> Optional[Tuple[int, float, float, float, float]]:
    """
    Extract key metrics from one simulation case:
    - iteration count
    - pressure outlet
    - head loss (p_in - p_out)
    - mass flow rate
    - max yPlus
    """
    try:
        files = get_required_files(case_dir)
        it      = read_value_from_last_line(files['solver'], 0)
        pin     = read_value_from_last_line(files['p_in'], -1)
        pout    = read_value_from_last_line(files['p_out'], -1)
        mrate   = read_value_from_last_line(files['m_flow'], -1)
        y_plus  = read_value_from_last_line(files['y_plus'], -2)
        return int(it), pout, (pin - pout), mrate, y_plus
    except Exception as e:
        print(f"[ERROR] Failed to extract data from {case_dir.name}")
        import traceback
        traceback.print_exc()
        return None

# ────────────────────────────────
# Plotting
# ────────────────────────────────

def plot_data(case_dir: Path, name: str) -> None:
    import matplotlib.pyplot as plt

    try:
        base = case_dir / "postProcessing"
        pout  = read_values(safe_glob(base, "averageTotalPressureOutlet1/0/surfaceFieldValue.dat"))
        mflow = read_values(safe_glob(base, "massFlowRateOutlet1/0/surfaceFieldValue.dat"))
        yplus = read_values(safe_glob(base, "yPlus1/0/yPlus.dat"), value_col=-2)
    except Exception:
        print(f"[ERROR] Failed to read plot data for {name}")
        return

    if any(len(x[0]) <= 1 for x in [pout, mflow, yplus]):
        print(f"[WARN] Not enough data to plot for {name}")
        return

    plt.figure(figsize=(10, 6))

    plt.subplot(3, 1, 1)
    plt.plot(*pout, label="Pressure Outlet", color="blue", linestyle='--', marker='o')
    plt.title(f"{name} - Pressure Outlet")
    plt.xlabel("Iterations")
    plt.ylabel("Pressure (Pa)")
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(*mflow, label="Mass Flow", color="green", linestyle='--', marker='o')
    plt.title(f"{name} - Mass Flow Rate")
    plt.xlabel("Iterations")
    plt.ylabel("kg/s")
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(*yplus, label="yPlus", color="red", linestyle='--', marker='o')
    plt.title(f"{name} - Max yPlus")
    plt.xlabel("Iterations")
    plt.ylabel("yPlus")
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# ────────────────────────────────
# Utility: Natural Sorting
# ────────────────────────────────

def natural_key(path: Path) -> List:
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', path.name)]

# ────────────────────────────────
# Geometry Parsing
# ────────────────────────────────

def parse_case_name(name: str) -> Tuple[Optional[str], Optional[str], Optional[Tuple[int, int]]]:
    """
    Identify case geometry and parameters from directory name.
    Returns (geometry, key, (Mb, Mw)) or (None, None, None) if unknown.
    """
    parts = name.split("-")
    try:
        if name.startswith("ELL-"):
            _, Mw, Mb, Kx, Ky, r, L, t = parts
            return "ELL", f"{Kx}_{Ky}_{r}_{L}", (int(Mb), int(Mw))
        elif name.startswith("IN-"):
            _, Mw, L, t = parts
            return "IN", L, (1, int(Mw))
        elif name.startswith("RAD-"):
            _, Mw, Mb, K, L, t = parts
            return "RAD", f"{K}_{L}", (int(Mb), int(Mw))
    except ValueError:
        pass
    return None, None, None

# ────────────────────────────────
# Main Workflow
# ────────────────────────────────

def main(write: bool = False, plot: bool = False) -> None:
    if not OUTPUTS_DIR.exists():
        print(f"[ERROR] Outputs directory not found: {OUTPUTS_DIR}")
        sys.exit(1)

    cases = sorted([p for p in OUTPUTS_DIR.iterdir() if p.is_dir()], key=natural_key)
    print(f"Found {len(cases)} cases.")

    results: Dict[Tuple[str, str], Dict[Tuple[int, int], Dict[str, float]]] = {}

    for case in cases:
        name = case.name
        print(f"\nProcessing: {name}")

        geom, key, dims = parse_case_name(name)
        if not geom:
            print(f"[WARN] Skipping unknown format: {name}")
            continue

        data = extract_case_data(case)
        if data is None:
            continue

        n_iter, pt_out, head_loss, massflow, max_yplus = data

        print(f"  - Iterations:      {n_iter}")
        print(f"  - Pressure Out:    {pt_out:.4f} Pa")
        print(f"  - Head Loss:       {head_loss:.4f} Pa")
        print(f"  - Mass Flow Rate:  {massflow:.4f} kg/s")
        print(f"  - Max yPlus:       {max_yplus:.2f}")
        if max_yplus > 5:
            print("  WARNING: yPlus > 5!")
        elif max_yplus >= 1:
            print("  WARNING: yPlus between 1 and 5.")

        results.setdefault((geom, key), {})[dims] = {
            "head_loss": head_loss,
            "massflow":  massflow,
            "max_yplus": max_yplus,
        }

        if plot:
            plot_data(case, name)

    if not write:
        print("\nDone. Use --write to export CSV files.")
        return

    POSTPROC_DIR.mkdir(exist_ok=True)
    import pandas as pd

    for (geom, key), data in results.items():
        Mb_vals = sorted(set(mb for mb, _ in data))
        Mw_vals = sorted(set(mw for _, mw in data))
        idx = pd.Index(Mb_vals, name="Mb")
        cols = pd.Index(Mw_vals, name="Mw")

        df_head = pd.DataFrame(index=idx, columns=cols, dtype=float)
        df_flow = df_head.copy()
        df_yplus = df_head.copy()

        for (mb, mw), vals in data.items():
            df_head.at[mb, mw] = vals["head_loss"]
            df_flow.at[mb, mw] = vals["massflow"]
            df_yplus.at[mb, mw] = vals["max_yplus"]

        base = f"{geom}_{key}"
        df_head.to_csv(POSTPROC_DIR / f"{base}_head_losses.csv")
        df_flow.to_csv(POSTPROC_DIR / f"{base}_massflow.csv")
        df_yplus.to_csv(POSTPROC_DIR / f"{base}_max_yplus.csv")
        print(f"[OK] Saved CSVs for {geom}/{key}")

    print("\nAll processing complete.")

# ────────────────────────────────
# Script Entry Point
# ────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    main(write=args.write, plot=args.plot)
