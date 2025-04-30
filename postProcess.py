#!/usr/bin/env python3
"""
Post-process OpenFOAM outputs: extract head losses, mass flow rates, and max yPlus
from each case under ./outputs/<caseName>/postProcessing, then assemble and
write CSV tables under ./postProcessingOutputs.
"""

import sys
from pathlib import Path
import glob

import numpy as np
import pandas as pd

# ─── Configuration ─────────────────────────────────────────────────────────────
OUTPUTS_DIR = Path("outputs")
POSTPROC_DIR = Path("postProcessingOutputs")
POSTPROC_DIR.mkdir(exist_ok=True)

# ─── I/O HELPERS ───────────────────────────────────────────────────────────────

def read_last_line(path: Path) -> str:
    """
    Return the last non-empty, non-comment line from `path`.
    """
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    # scan backwards in file
    with path.open() as f:
        f.seek(0, 2)
        pos = f.tell()
        buffer = ""
        while pos > 0:
            pos -= 1
            f.seek(pos)
            char = f.read(1)
            if char == "\n":
                line = buffer[::-1].strip()
                if line and not line.startswith("#"):
                    return line
                buffer = ""
            else:
                buffer += char
        # fallback: read all lines
        f.seek(0)
        valid = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
        if valid:
            return valid[-1]
    raise ValueError(f"No data in {path}")

def read_value(path: Path, index: int = -1) -> float:
    """
    Read a float from the last valid line of `path`, taking token at `index`.
    """
    line = read_last_line(path)
    tokens = line.split()
    if len(tokens) < abs(index) + 1:
        raise ValueError(f"Not enough tokens in {path}: {tokens}")
    return float(tokens[index])

# ─── EXTRACTION ────────────────────────────────────────────────────────────────

def extract_case_data(case_dir: Path):
    """
    From one case directory, extract:
      - iteration (first token of solverInfo.dat)
      - head loss (p_inlet - p_outlet)
      - mass flow (last token of massFlowRateOutlet1)
      - max yPlus (second-last token of yPlus.dat)
    Returns (iteration, head_loss, massflow, max_yplus) or None on error.
    """

    try:
        base = case_dir / "postProcessing"
        # find each file
        solver = next(base.glob("solverInfo1/0/solverInfo.dat"))
        p_in   = next(base.glob("averageTotalPressureInlet1/0/surfaceFieldValue.dat"))
        p_out  = next(base.glob("averageTotalPressureOutlet1/0/surfaceFieldValue.dat"))
        mflow  = next(base.glob("massFlowRateOutlet1/0/surfaceFieldValue.dat"))
        yplus  = next(base.glob("yPlus1/0/yPlus.dat"))
        # read data
        it      = read_value(solver, 0)
        pin     = read_value(p_in,   -1)
        pout    = read_value(p_out,  -1)
        mrate   = read_value(mflow,  -1)
        y_plus  = read_value(yplus,  -2)
        return it, (pin - pout), mrate, y_plus

    except Exception as e:
        import traceback
        print(f"[ERROR] Failed to process {case_dir.name}")
        traceback.print_exc()
        return None


# ─── MAIN WORKFLOW ────────────────────────────────────────────────────────────

def main():
    if not OUTPUTS_DIR.is_dir():
        print(f"[ERROR] Outputs folder not found: {OUTPUTS_DIR}")
        sys.exit(1)

    # gather all case subdirectories
    cases = [p for p in OUTPUTS_DIR.iterdir() if p.is_dir()]
    print(f"Found {len(cases)} cases in {OUTPUTS_DIR}")

    # data structure: { (geom, params_key): { (Mb, Mw): (head, mflow, yplus) } }
    results = {}

    for case in sorted(cases):
        name = case.name
        print(f"\nProcessing case: {name}")
        # determine geometry type & parameters key
        if name.startswith("ELL-"):
            _, Mw, Mb, Kx, Ky, r, L, t = name.split("-")
            geom, key, dims = "ELL", f"{Kx}_{Ky}_{r}_{L}", (int(Mb), int(Mw))
        elif name.startswith("IN-"):
            _, Mw, L, t = name.split("-")
            geom, key, dims = "IN", L, (1, int(Mw))
        elif name.startswith("RAD-"):
            _, Mw, Mb, K, L, t = name.split("-")
            geom, key, dims = "RAD", f"{K}_{L}", (int(Mb), int(Mw))
        else:
            print(f"[WARN] Unknown pattern: {name}, skipping")
            continue
        
        data = extract_case_data(case)
        if data is None:
            continue

        
        n_iter, head_loss, massflow, max_yplus = data

        # Reporting
        print(f"  - Iterations              : {int(n_iter)}")
        print(f"  - Mass flow rate (kg/s)   : {massflow:.4f}")
        print(f"  - Head loss (Pa)          : {head_loss:.2f}")
        print(f"  - Max yPlus               : {max_yplus:.2f}")
        if max_yplus > 5.0:
            print("  WARNING: yPlus > 5 detected!")
        elif max_yplus >= 1.0:
            print("  WARNING: yPlus not fully below 1!")
        
        results.setdefault((geom, key), {})[dims] = {
            "head_loss": head_loss,
            "massflow":  massflow,
            "max_yplus": max_yplus
        }


    # for each geometry+params, build and save DataFrames
    for (geom, key), table in results.items():
        # build index/columns
        Mb_vals = sorted({mb for mb, _ in table})
        Mw_vals = sorted({mw for _, mw in table})
        idx = pd.Index(Mb_vals, name="Mb")
        cols = pd.Index(Mw_vals, name="Mw")

        # initialize empty DataFrames
        df_head  = pd.DataFrame(index=idx, columns=cols, dtype=float)
        df_mflow = df_head.copy()
        df_yplus = df_head.copy()

        # fill
        for (mb, mw), vals in table.items():
            df_head.at[mb, mw]  = vals["head_loss"]
            df_mflow.at[mb, mw] = vals["massflow"]
            df_yplus.at[mb, mw] = vals["max_yplus"]

        # output CSVs
        base = f"{geom}_{key}"
        df_head.to_csv(POSTPROC_DIR / f"{base}_head_losses.csv")
        df_mflow.to_csv(POSTPROC_DIR / f"{base}_massflow.csv")
        df_yplus.to_csv(POSTPROC_DIR / f"{base}_max_yplus.csv")

        print(f"[OK] Wrote {geom}/{key} tables to {POSTPROC_DIR}")

    print("\nAll post-processing complete.")

if __name__ == "__main__":
    main()
