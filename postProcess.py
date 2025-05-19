#!/usr/bin/env python3
"""
Main Features:
- Parses OpenFOAM post-processing output files for multiple simulation cases.
- Extracts key metrics: iteration count, outlet pressure, head loss, mass flow rate, and maximum yPlus.
- Supports natural sorting of case directories.
- Handles multiple geometry types (ELL, IN, RAD) based on directory naming conventions.
- Writes results to CSV tables for head loss, mass flow, and yPlus, organized by geometry and case parameters.
- Optionally generates and saves plots for each case, including relative error curves.
- Provides detailed logging and error handling for missing or malformed files.
Usage:
    python postProcess.py [--write] [--plot]
Arguments:
    --write   Write results to CSV files.
    --plot    Display plots for each case.
Dependencies:
    - Python 3.6+
    - matplotlib
    - pandas (only if --write is used)
"""

import sys
import argparse
import re
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Union
import matplotlib.pyplot as plt
import logging
import traceback
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
    Reads the last non-empty, non-comment line from a file.

    Args:
        path (Path): The path to the file to read.

    Returns:
        str: The last valid line (not empty and not starting with '#') from the file.

    Raises:
        ValueError: If no valid data line is found in the file.
    """
    with path.open("r") as f:
        for line in reversed(f.readlines()):
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                return stripped
    raise ValueError(f"No valid data in {path}")

def read_value_from_last_line(path: Path, index: int = -1) -> float:
    """
    Reads a floating-point value from a specific token in the last line of a file.

    Args:
        path (Path): The path to the file from which to read.
        index (int, optional): The index of the token to extract from the last line.
            Defaults to -1 (the last token).

    Returns:
        float: The value of the specified token converted to a float.

    Raises:
        ValueError: If the last line does not contain enough tokens to satisfy the index.
    """
    line = read_last_line(path)
    tokens = line.split()
    if len(tokens) < abs(index) + 1:
        raise ValueError(f"Not enough tokens in {path}: {tokens}")
    return float(tokens[index])

def read_values(file_path: Path, iter_col: int = 0, value_col: int = -1) -> Tuple[List[float], List[float]]:
    """
    Reads numerical data from a file, extracting iteration and value columns.

    Args:
        file_path (Path): Path to the file containing data.
        iter_col (int, optional): Index of the column to use as the iteration or x-axis values. Defaults to 0.
        value_col (int, optional): Index of the column to use as the value or y-axis values. Defaults to -1 (last column).

    Returns:
        Tuple[List[float], List[float]]: Two lists containing the iteration values and the corresponding data values.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If a line in the file is malformed or cannot be converted to float.

    Notes:
        - Lines starting with '#' or empty lines are ignored.
        - Each valid line is split by whitespace.
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
    Searches for the first file matching the given pattern in the specified base directory.

    Args:
        base (Path): The base directory to search in.
        pattern (str): The glob pattern to match files.

    Returns:
        Path: The first Path object matching the pattern.

    Raises:
        FileNotFoundError: If no file matching the pattern is found in the base directory.
    """
    try:
        return next(base.glob(pattern))
    except StopIteration:
        raise FileNotFoundError(f"No file matching {pattern} in {base}")

def get_required_files(case_dir: Path) -> Dict[str, Path]:
    """
    Collects and returns the required post-processing file paths for a given OpenFOAM case directory.

    Args:
        case_dir (Path): The root directory of the OpenFOAM case.

    Returns:
        Dict[str, Path]: A dictionary mapping descriptive keys to the corresponding file paths for:
            - 'solver': Solver information data file.
            - 'p_in': Average total pressure at the inlet.
            - 'p_out': Average total pressure at the outlet.
            - 'm_flow': Mass flow rate at the outlet.
            - 'y_plus': yPlus values.

    Note:
        Uses the `safe_glob` function to locate files within the 'postProcessing' subdirectory.
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
    Extracts simulation data from a given case directory.

    This function retrieves specific values from required files within the provided
    case directory. It reads the last line (or specified line) from each file to extract:
    - The iteration number from the solver file.
    - The inlet pressure from the p_in file.
    - The outlet pressure from the p_out file.
    - The mass flow rate from the m_flow file.
    - The y-plus value from the y_plus file.

    Returns a tuple containing:
        (iteration, outlet_pressure, pressure_difference, mass_flow_rate, y_plus)

    If any error occurs during extraction, logs the error and returns None.

    Args:
        case_dir (Path): The path to the case directory containing the required files.

    Returns:
        Optional[Tuple[int, float, float, float, float]]: 
            A tuple with (iteration, outlet_pressure, pressure_difference, mass_flow_rate, y_plus),
            or None if extraction fails.
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
        print("")
        logger.error(f"Failed to extract data from {case_dir.name}")
        traceback.print_exc()
        return None

# ────────────────────────────────
# Plotting
# ────────────────────────────────

def plot_data(case_dir: Path, name: str, plot:bool=False, write:bool=False, postproc_dir:Path=None) -> None:
    """
    Plots and/or saves diagnostic data from OpenFOAM post-processing results.
    This function reads simulation output data (total pressure, mass flow rate, and yPlus)
    from specified directories, generates plots for each quantity along with their relative errors,
    and optionally displays or saves the plots.

    Args:
        case_dir (Path): Path to the OpenFOAM case directory.
        name (str): Name identifier for the plot (used in titles and filenames).
        plot (bool, optional): If True, displays the plot interactively. Defaults to False.
        write (bool, optional): If True, saves the plot as a PNG file in `postproc_dir`. Defaults to False.
        postproc_dir (Path, optional): Directory where the plot image will be saved if `write` is True. Defaults to None.
    
    Raises:
        None explicitly. Logs errors and warnings if data cannot be read or is insufficient.
    
    Notes:
        - Requires matplotlib for plotting.
        - Assumes helper functions `read_values`, `safe_glob`, and a `logger` are defined elsewhere.
        - Plots three subplots: Outlet Total Pressure, Outlet Mass Flow Rate, and Max yPlus,
          each with their relative error on a secondary y-axis.
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

    def rel_error(values, eps=1e-10) -> List[float]:
        """
        Calculate the relative percentage error between consecutive values in a list.

        Args:
            values (List[float]): A list of numerical values to compute relative errors for.
            eps (float, optional): A small value to prepend to the result list. Defaults to 1e-10.

        Returns:
            List[float]: A list where the first element is `eps`, followed by the relative percentage errors
            between each consecutive pair of values in `values`. If the previous value is zero, the error is set to 0.
        """

        return [eps] + [abs((v2 - v1) / v1)*100 if v1 != 0 else 0 for v1, v2 in zip(values[:-1], values[1:])]
    
    _, axs = plt.subplots(3, 1, figsize=(10, 8))

    data_pairs = [
        (pout, "Outlet Total Pressure", "Pressure (Pa/m)", "blue"),
        (mflow, "Outlet Mass Flow Rate", "Mass Flow Rate (kg/s/m)", "green"),
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
        # ax2.set_yscale('log')
        # ax2.set_ylim(-0.01, 0.5)
        ax2.tick_params(axis='y', labelcolor="gray")

    plt.tight_layout()

    if write:
        postproc_dir.mkdir(parents=True, exist_ok=True)
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
    Generate a key for natural sorting of file paths.

    This function splits the file name into a list of integers and lowercase strings,
    allowing for natural sorting (e.g., "file2" comes before "file10").

    Args:
        path (Path): The file path whose name will be used for sorting.

    Returns:
        List: A list of integers and strings representing the split components of the file name,
              suitable for use as a sorting key.
    """
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', path.name)]

# ────────────────────────────────
# Geometry Parsing
# ────────────────────────────────

def parse_case_name(name: str) -> Tuple[Optional[str], Optional[str], Optional[Tuple[int, int]]]:
    """
    Parses a case name string and extracts case type, parameters, and a tuple of integers.

    Args:
        name (str): The case name string to parse. Expected formats:
            - "ELL-<Mw>-<Mb>-<Kx>-<Ky>-<t>-<r>"
            - "IN-<Mw>"
            - "RAD-<Mw>-<Mb>-<K>-<t>"
            - "AP-<Mw>-<Mb>-<Kx>-<Ky>-<t>-<r>"

    Returns:
        Tuple[Optional[str], Optional[str], Optional[Tuple[int, int]]]:
            - A tuple containing:
                - The case type as a string ("ELL", "IN", or "RAD"), or None if parsing fails.
                - A string of extracted parameters, or an empty string/None if not applicable.
                - A tuple of two integers (Mb, Mw) or (0, Mw), or None if parsing fails.

    Notes:
        - Returns (None, None, None) if the input does not match any expected format or if parsing fails.
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
        elif name.startswith("AP-"):
            _, Mw, Mb, Kx, Ky, t, r = parts
            return "AP", f"{Kx}-{Ky}-{t}-{r}", (int(Mb), int(Mw))
    except ValueError:
        pass
    return None, None, None

# ────────────────────────────────
# Process Case
# ────────────────────────────────

def process_single_case(case: Path, postproc_dir: Path, plot: bool, write: bool) -> Optional[Tuple[Tuple[str, str], Tuple[int, int], Dict[str, Union[float, str]]]]:
    """
    Processes a single simulation case directory, extracting relevant data and optionally plotting or writing results.
    Args:
        case (Path): Path to the simulation case directory.
        postproc_dir (Path): Path to the directory where post-processing outputs should be saved.
        plot (bool): If True, generate plots for the case.
        write (bool): If True, write post-processed data to files.
    Returns:
        Optional[Tuple[Tuple[str, str], Tuple[int, int], Dict[str, Union[float, str]]]]:
            - ((geometry, key), (Mb, Mw), out_data) if the case is successfully processed, where:
                - geometry (str): Geometry identifier parsed from the case name.
                - key (str): Additional key parsed from the case name.
                - Mb (int): First dimension parsed from the case name.
                - Mw (int): Second dimension parsed from the case name.
                - out_data (dict): Dictionary containing extracted and computed case data, including:
                    - "geometry": Geometry identifier.
                    - "key": Key identifier.
                    - "Mw": Second dimension.
                    - "Mb": First dimension.
                    - "pt_out": Outlet total pressure.
                    - "head_loss": Head loss value.
                    - "massflow": Mass flow rate.
                    - "max_yplus": Maximum y+ value.
                    - "iterations": Number of solver iterations.
                    - "case_name": Name of the case directory.
            - None if the case format is unknown or data extraction fails.
    Notes:
        - Relies on external functions: parse_case_name, extract_case_data, plot_data.
        - Logs a warning and skips cases with unknown naming formats.
    """
    name = case.name

    geom, key, dims = parse_case_name(name)

    if not geom:
        print("")
        logger.warning(f"Skipping unknown format: {name}")
        return None
    data = extract_case_data(case)
    if data is None:
        return None

    n_iter, pt_out, head_loss, massflow, max_yplus = data
    
    out_data = {
        "geometry": geom,
        "key": key,
        "Mw": dims[1],
        "Mb": dims[0],
        "pt_out":     pt_out,
        "head_loss":  head_loss,
        "massflow":   massflow,
        "max_yplus":  max_yplus,
        "iterations": n_iter,
        "case_name":  name
    }

    if plot or write:
        plot_data(case, name, plot, write, postproc_dir)

    return (geom, key), dims, out_data
    
# ────────────────────────────────
# Main Workflow
# ────────────────────────────────

def main(write: bool = False, plot: bool = False) -> None:
    """
    Main entry point for post-processing simulation results.
    This function processes simulation output directories, aggregates results, and optionally writes summary CSV files and generates plots.
    Args:
        write (bool, optional): If True, writes processed results and summaries to CSV files. Defaults to False.
        plot (bool, optional): If True, generates plots for each case during processing. Defaults to False.
    Workflow:
        - Checks for the existence of the base outputs directory.
        - Iterates over all subdirectories (simulation groups) in the outputs directory.
        - For each group, processes all simulation cases, extracting relevant metrics.
        - Logs key statistics for each case, including head loss, mass flow rate, and max yPlus.
        - Aggregates results by geometry and key, organizing them by mesh parameters.
        - If `write` is True, saves head loss, mass flow, and max yPlus data as CSV files for each group.
        - Also writes a summary CSV of all cases if `write` is True.
        - Logs progress and warnings (e.g., if yPlus exceeds thresholds).
        - Reports total runtime upon completion.
    Raises:
        SystemExit: If the base outputs directory or required subdirectories are missing.
    """
    base_outputs = Path("outputs")
    if not base_outputs.exists():
        logger.error(f"Base outputs directory not found: {base_outputs}")
        sys.exit(1)

    all_subdirs = [d for d in base_outputs.iterdir() if d.is_dir()]
    if not all_subdirs:
        logger.error("No subdirectories found in outputs/.")
        sys.exit(1)

    for subdir in all_subdirs:
        summary_data = []
        output_dir = subdir
        postproc_dir = Path("postProcessingOutputs") / subdir.name
        
        cases = sorted([p for p in output_dir.iterdir() if p.is_dir()], key=natural_key)
        if not cases:
            print("")
            logger.info(f"No cases found in {output_dir}")
            continue
        print("")
        logger.info(f">>> Processing '{subdir.name}' with {len(cases)} cases")

        results: Dict[Tuple[str, str], Dict[Tuple[int, int], Dict[str, float]]] = {}

        for case in cases:

            single_case_output = process_single_case(case, postproc_dir, plot, write)
            
            if not single_case_output:
                    continue
            (geom, key), dims, data = single_case_output

            print("")
            logger.info(f"Processing: {data['case_name']}")
            logger.info(f"   - Iterations:     {data['iterations']}")
            logger.info(f"   - Pressure Out:   {data['pt_out']:.4f} Pa/m")
            logger.info(f"   - Head Loss:      {data['head_loss']:.4f} Pa/m")
            logger.info(f"   - Mass Flow Rate: {data['massflow']:.4f} kg/s/m")
            logger.info(f"   - Max yPlus:      {data['max_yplus']:.2f}")
            if data['max_yplus'] > 5:
                logger.warning("yPlus > 5!")
            elif data['max_yplus'] >= 1:
                logger.warning("yPlus between 1 and 5.")

            results.setdefault((geom, key), {})[dims] = {
                "head_loss": data['head_loss'],
                "massflow":  data['massflow'],
                "max_yplus": data['max_yplus'],
            }
            summary_data.append(data)

        if write and results:
            import pandas as pd
            postproc_dir.mkdir(parents=True, exist_ok=True)

            print("")
            for (geom, key), data_dict in results.items():
                Mb_vals = sorted(set(mb for mb, _ in data_dict))
                Mw_vals = sorted(set(mw for _, mw in data_dict))
                idx = pd.Index(Mb_vals, name="Mb")
                cols = pd.Index(Mw_vals, name="Mw")

                df_head = pd.DataFrame(index=idx, columns=cols, dtype=float)
                df_flow = df_head.copy()
                df_yplus = df_head.copy()

                for (mb, mw), vals in data_dict.items():
                    df_head.at[mb, mw] = vals["head_loss"]
                    df_flow.at[mb, mw] = vals["massflow"]
                    df_yplus.at[mb, mw] = vals["max_yplus"]

                base = f"{geom}" if geom == "IN" else f"{geom}-{key}"
                df_head.to_csv(postproc_dir / f"{base}_head_losses.csv")
                df_flow.to_csv(postproc_dir / f"{base}_massflow.csv")
                df_yplus.to_csv(postproc_dir / f"{base}_max_yplus.csv")
                logger.info(f"Saved CSVs for {base} in {postproc_dir}")       
        
        if write and summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(postproc_dir/ "all_cases_summary.csv", index=False)
            logger.info(f"Saved all-cases summary to {postproc_dir}/all_cases_summary.csv")
    
    print("")
    logger.info("All processing complete.")
    logger.info(f"Total runtime: {time.time() - start_time:.2f} seconds")

# ────────────────────────────────
# Script Entry Point
# ────────────────────────────────

if __name__ == "__main__":
    args = parse_args()
    main(write=args.write, plot=args.plot)
