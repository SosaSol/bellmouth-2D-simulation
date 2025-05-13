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
import matplotlib.pyplot as plt
import logging
import time

logging.basicConfig(level=logging.INFO, format='[%(levelname)s]\t%(message)s')
logger = logging.getLogger(__name__)
logging.addLevelName(logging.WARNING, "WARN")   # Override the default 'WARNING' level name

start_time = time.time()

# ────────────────────────────────
# Argument Parsing
# ────────────────────────────────

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
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

    Args:
        path (Path): Path to the file.
    Returns:
        str: The last valid line.
    Raises:
        ValueError: If no valid line is found.
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

    Args:
        path (Path): Path to the file.
        index (int): Index of the token to extract (0-based, negative for reverse).
    Returns:
        float: The extracted value.
    Raises:
        ValueError: If the index is out of range or if the file is empty.
    """
    line = read_last_line(path)
    tokens = line.split()
    if len(tokens) < abs(index) + 1:
        raise ValueError(f"Not enough tokens in {path}: {tokens}")
    return float(tokens[index])

def read_values(file_path: Path, iter_col: int = 0, value_col: int = -1) -> Tuple[List[float], List[float]]:
    """
    Read columns of iteration and value data from a file.
    Ignores empty lines and comments.

    Args:
        file_path (Path): Path to the file.
        iter_col (int): Column index for iterations (0-based).
        value_col (int): Column index for values (0-based, negative for reverse).
    Returns:
        Tuple[List[float], List[float]]: Two lists containing iterations and values.
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is malformed or if the indices are out of range.
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
    Safely glob for a file, raising an error if not found.
    Args:
        base (Path): Base directory to search in.
        pattern (str): Glob pattern to match.
    Returns:
        Path: The first matching file.
    Raises:
        FileNotFoundError: If no matching file is found.
    """
    try:
        return next(base.glob(pattern))
    except StopIteration:
        raise FileNotFoundError(f"No file matching {pattern} in {base}")

def get_required_files(case_dir: Path) -> Dict[str, Path]:
    """
    Collect paths to required post-processing output files.
    Args:
        case_dir (Path): Directory containing the case data.
    Returns:
        dict (Dict[str, Path]): Dictionary mapping file identifiers to paths.
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
    Args:
        case_dir (Path): Directory containing the case data.
    Returns:
        case_data (Optional[Tuple[int, float, float, float, float]]): Tuple containing:
            - iteration count
            - pressure outlet
            - head loss
            - mass flow rate
            - max yPlus
    Raises:
        Exception: If any required file is missing or malformed.
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
        logger.error(f"Failed to extract data from {case_dir.name}")
        import traceback
        traceback.print_exc()
        return None

# ────────────────────────────────
# Plotting
# ────────────────────────────────

def plot_data(case_dir: Path, name: str, plot:bool=False, write:bool=False, postproc_dir:Path=None) -> None:
    """
    Plot post-processing data for a given case.
    Generates plots for pressure outlet, mass flow rate, and yPlus.
    Saves plots to the case directory if write is True.
    
    Args:
        case_dir (Path): Directory containing the case data.
        name (str): Name of the case.
        plot (bool): Whether to display the plots.
        write (bool): Whether to save the plots as PNG files.
        postproc_dir (Path): Directory to save the plots.
    """

    try:
        base = case_dir / "postProcessing"
        pout  = read_values(safe_glob(base, "averageTotalPressureOutlet1/0/surfaceFieldValue.dat"))
        mflow = read_values(safe_glob(base, "massFlowRateOutlet1/0/surfaceFieldValue.dat"))
        yplus = read_values(safe_glob(base, "yPlus1/0/yPlus.dat"), value_col=-2)
    except Exception:
        logger.error(f"Failed to read plot data for {name}")
        return

    if any(len(x[0]) <= 1 for x in [pout, mflow, yplus]):
        logger.warning(f"Not enough data to plot for {name}")
        return

    def rel_error(values) -> List[float]:
        """Compute relative error in % between consecutive steps.
        Args:
            values (list): List of values.
        Returns:
            list: List of relative errors in %.
        """
        return [0] + [abs((v2 - v1) / v1)*100 if v1 != 0 else 0 for v1, v2 in zip(values[:-1], values[1:])]
    
    _, axs = plt.subplots(3, 1, figsize=(10, 8))

    data_pairs = [
        (pout, "Total Pressure Outlet", "Pressure (Pa)", "blue"),
        (mflow, "Outlet Mass Flow Rate", "Outlet Mass Flow Rate (kg/s)", "green"),
        (yplus, "Max yPlus", "Max y+", "red")
    ]

    for ax, (data, title, ylabel, color) in zip(axs, data_pairs):
        x, y = data
        ax.plot(x, y, label=title, color=color, linestyle='--', marker='o')
        ax.set_title(f"{name} - {title}")
        ax.set_xlabel("Iterations")
        ax.set_ylabel(ylabel)
        ax.grid(True)
        # remove offset from y-axis
        ax.ticklabel_format(style='plain', axis='y', useOffset=False)

        # Plot relative error on second y-axis
        ax2 = ax.twinx()
        rel_err = rel_error(y)
        ax2.plot(x[2:], rel_err[2:], color="gray", linestyle=":", marker='o', fillstyle='none' , label="Rel. Error")
        ax2.set_ylabel("Rel. Error (%)", color="gray")
        ax2.tick_params(axis='y', labelcolor="gray")

    plt.tight_layout()

    if write:
        plt.savefig(postproc_dir / f"{name}.png")
        plt.close()
        logger.info(f"Saved plot for {name} to {postproc_dir}")
    if plot:
        plt.show()
    else:
        plt.close()
# ────────────────────────────────
# Utility: Natural Sorting
# ────────────────────────────────

def natural_key(path: Path) -> List:
    """
    Generate a natural sorting key for a path.
    Args:
        path (Path): Path to the file or directory.
    Returns:
        list (List): A list of components for natural sorting.
    """
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', path.name)]

# ────────────────────────────────
# Geometry Parsing
# ────────────────────────────────

def parse_case_name(name: str) -> Tuple[Optional[str], Optional[str], Optional[Tuple[int, int]]]:
    """
    Identify case geometry and parameters from directory name.
    Returns (geometry, key, (Mb, Mw)) or (None, None, None) if unknown.
    Args:
        name (str): Name of the case directory.
    Returns:
        tuple (Tuple[Optional[str], Optional[str], Optional[Tuple[int, int]]]): Tuple containing:
            - geometry type (e.g., "ELL", "IN", "RAD")
            - key (e.g., "Kx-Ky-t-r")
            - dimensions (Mb, Mw)
    Raises:
        ValueError: If the name format is unrecognized.
    """
    parts = name.split("-")
    try:
        if name.startswith("ELL-"):
            _, Mw, Mb, Kx, Ky, t, r = parts
            return "ELL", f"{Kx}-{Ky}-{t}-{r}", (int(Mb), int(Mw))
        elif name.startswith("IN-"):
            _, Mw, = parts
            return "IN","", (0, int(Mw))
        elif name.startswith("RAD-"):
            _, Mw, Mb, K, t = parts
            return "RAD", f"{K}-{t}", (int(Mb), int(Mw))
    except ValueError:
        pass
    return None, None, None

# ────────────────────────────────
# Main Workflow
# ────────────────────────────────

def main(write: bool = False, plot: bool = False) -> None:
    """
    Main function to process OpenFOAM outputs.
    Args:
        write (bool): Whether to write results to CSV files.
        plot (bool): Whether to generate plots.
    """
    base_outputs = Path("outputs")
    if not base_outputs.exists():
        logger.error(f"Base outputs directory not found: {base_outputs}")
        sys.exit(1)

    all_subdirs = [d for d in base_outputs.iterdir() if d.is_dir()]
    if not all_subdirs:
        logger.error("No subdirectories found in outputs/.")
        sys.exit(1)

    summary_data = []

    for subdir in all_subdirs:
        output_dir = subdir
        postproc_dir = Path("postProcessingOutputs") / subdir.name
        postproc_dir.mkdir(parents=True, exist_ok=True)

        cases = sorted([p for p in output_dir.iterdir() if p.is_dir()], key=natural_key)
        print("")
        logger.info(f">>> Processing '{subdir.name}' with {len(cases)} cases")

        results: Dict[Tuple[str, str], Dict[Tuple[int, int], Dict[str, float]]] = {}

        for case in cases:
            name = case.name
            print("")
            logger.info(f"Processing: {name}")

            geom, key, dims = parse_case_name(name)
            if not geom:
                logger.warning(f"Skipping unknown format: {name}")
                continue

            data = extract_case_data(case)
            if data is None:
                continue

            n_iter, pt_out, head_loss, massflow, max_yplus = data

            logger.info(f"  - Iterations:      {n_iter}")
            logger.info(f"  - Pressure Out:    {pt_out:.4f} Pa")
            logger.info(f"  - Head Loss:       {head_loss:.4f} Pa")
            logger.info(f"  - Mass Flow Rate:  {massflow:.4f} kg/s")
            logger.info(f"  - Max yPlus:       {max_yplus:.2f}")
            if max_yplus > 5:
                logger.warning("yPlus > 5!")
            elif max_yplus >= 1:
                logger.warning("yPlus between 1 and 5.")

            results.setdefault((geom, key), {})[dims] = {
                "head_loss": head_loss,
                "massflow":  massflow,
                "max_yplus": max_yplus,
            }

            summary_data.append({
                "geometry": geom,
                "key": key,
                "Mb": dims[0],
                "Mw": dims[1],
                "head_loss": head_loss,
                "massflow": massflow,
                "max_yplus": max_yplus,
                "case_name": name
            })

            if plot or write:
                plot_data(case, name, plot, write, postproc_dir)

        if write:
            import pandas as pd

            print("")
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

                if geom == "IN":
                    base = f"{geom}"
                else:
                    base = f"{geom}-{key}"
                
                df_head.to_csv(postproc_dir / f"{base}_head_losses.csv")
                df_flow.to_csv(postproc_dir / f"{base}_massflow.csv")
                df_yplus.to_csv(postproc_dir / f"{base}_max_yplus.csv")
                logger.info(f"Saved CSVs for {base} in {postproc_dir}")
            
        if write and summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(Path("postProcessingOutputs") / "all_cases_summary.csv", index=False)
            logger.info("Saved all-cases summary to postProcessingOutputs/all_cases_summary.csv")


    print("")
    logger.info("All processing complete.")
    logger.info(f"Total runtime: {time.time() - start_time:.2f} seconds")

# ────────────────────────────────
# Script Entry Point
# ────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    main(write=args.write, plot=args.plot)
