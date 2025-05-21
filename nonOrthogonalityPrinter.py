#!/usr/bin/env python3
"""
nonOrthogonalityPrinter.py
This script scans through a specified directory structure containing simulation output folders,
searches for "log.checkMesh" files, and extracts mesh non-orthogonality statistics (maximum and average values)
from each log file. It prints a color-coded summary of these values for each case and highlights the cases
with the highest maximum and average non-orthogonality.
Functions:
-----------
find_non_orthogonality(log_file_path: str) -> tuple:
    Parses the given log file to extract the maximum and average mesh non-orthogonality values.
    Returns a tuple (max_val, avg_val) as strings.
    Raises ValueError if the values cannot be found.
main():
    Main entry point of the script.
    - Accepts an optional command-line argument specifying the outputs directory (defaults to './outputs').
    - Iterates through each folder and subfolder, searching for "log.checkMesh" files.
    - Extracts and prints non-orthogonality statistics, color-coding the output based on severity.
    - Identifies and prints the cases with the highest maximum and average non-orthogonality.
Usage:
------
    python nonOrthogonalityPrinter.py [outputs_directory]
Dependencies:
-------------
    - natsort (optional, for natural sorting of folder names)
    - Standard Python libraries: os, logging, re, sys
Color Coding:
-------------
    - Red:   Max non-orthogonality >= 80
    - Yellow: Max non-orthogonality >= 70 and < 80
    - Default: Max non-orthogonality < 70
"""

import os
import logging
import re
import sys

try:
    from natsort import natsorted
except ImportError:
    def natsorted(seq):
        """
        Sorts a sequence of strings in natural (human) order.

        This function sorts strings containing numbers in a way that is intuitive to humans,
        so that, for example, 'item2' comes before 'item10'.

        Args:
            seq (iterable): A sequence of strings to be sorted.

        Returns:
            list: A new list containing the sorted strings in natural order.

        Example:
            >>> natsorted(['item10', 'item2', 'item1'])
            ['item1', 'item2', 'item10']
        """
        import re
        def atoi(text):
            return int(text) if text.isdigit() else text
        def natural_keys(text):
            return [atoi(c) for c in re.split(r'(\d+)', text)]
        return sorted(seq, key=natural_keys)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Add file handler to save logs to 'nonOrthogonalityPrinter.log'
file_handler = logging.FileHandler('nonOrthogonalityPrinter.log', mode='w')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(logging.Formatter('%(message)s'))
logging.getLogger().addHandler(file_handler)

def find_non_orthogonality(log_file_path: str) -> tuple:
    """
    Extracts the maximum and average mesh non-orthogonality values from an OpenFOAM log file.

    Args:
        log_file_path (str): The path to the OpenFOAM log file to be parsed.

    Returns:
        tuple: A tuple containing two strings:
            - max_val (str): The maximum mesh non-orthogonality value found in the log file.
            - avg_val (str): The average mesh non-orthogonality value found in the log file.

    Raises:
        ValueError: If the non-orthogonality values cannot be found in the log file.
        FileNotFoundError: If the specified log file does not exist.
        IOError: If there is an error reading the log file.

    Logs:
        An error message is logged if the file cannot be read.
    """
    pattern = re.compile(r"Mesh non-orthogonality Max:\s*([\d\.]+)\s*average:\s*([\d\.]+)")
    try:
        with open(log_file_path, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    max_val = match.group(1)
                    avg_val = match.group(2)
                    return max_val, avg_val
    except (FileNotFoundError, IOError) as e:
        logging.error(f"Error reading log file {log_file_path}: {e}")
    raise ValueError(f"Failed to find non-orthogonality values in the log file: {log_file_path}")

def main():
    """
    Main function to extract and report non-orthogonality values from OpenFOAM simulation output directories.
    This function processes a specified output directory (or a default 'outputs' directory if none is provided),
    searching for 'log.checkMesh' files within subdirectories. It extracts the maximum and average non-orthogonality
    values from these log files using the `find_non_orthogonality` function, and prints color-coded summaries to the console.
    Additionally, it logs the cases with the highest maximum and average non-orthogonality values.
    Command-line Arguments:
        sys.argv[1] (optional): Path to the outputs directory. If not provided, defaults to './outputs'.
    Logging:
        - Logs errors if the outputs directory does not exist.
        - Logs progress and summary information for each processed folder and subfolder.
        - Logs warnings if a log file cannot be processed.
    Console Output:
        - Prints color-coded non-orthogonality values for each case.
        - Highlights cases with the highest maximum and average values.
    Dependencies:
        - Requires the `find_non_orthogonality` function to extract values from log files.
        - Uses `natsorted` for natural sorting of directory names.
        - Uses `logging` for logging messages.
    """
    if len(sys.argv) > 1:
        outputs_dir = os.path.abspath(sys.argv[1])
    else:
        outputs_dir = os.path.join(os.getcwd(), "outputs")
    if not os.path.isdir(outputs_dir):
        logging.error(f"'{outputs_dir}' directory not found")
        return

    logging.info(f"Extracting non-orthogonality values in {outputs_dir}")

    for folder in natsorted(os.listdir(outputs_dir)):
        folder_path = os.path.join(outputs_dir, folder)
        if os.path.isdir(folder_path):
            print("")
            logging.info(f">> Processing {folder} folder")
            folder_results = []
            for subfolder in natsorted(os.listdir(folder_path)):
                subfolder_path = os.path.join(folder_path, subfolder)
                if os.path.isdir(subfolder_path):
                    log_file = os.path.join(subfolder_path, "log.checkMesh")
                    if os.path.isfile(log_file):
                        try:
                            max_val, avg_val = find_non_orthogonality(log_file)
                            if max_val and avg_val:
                                max_val_f = float(max_val)
                                avg_val_f = float(avg_val)
                                folder_results.append((subfolder, max_val_f, avg_val_f))
                                # Color output for terminal only
                                if max_val_f >= 80:
                                    color = "\033[91m"  # Red
                                elif max_val_f >= 70:
                                    color = "\033[93m"  # Orange/Yellow
                                else:
                                    color = "\033[0m"   # Default
                                reset = "\033[0m"
                                logging.info(f"{color}{subfolder:<20}\tMax={max_val_f:.1f}, Avg={avg_val_f:.1f}{reset}")
                        except Exception as e:
                            logging.warning(f"Could not process {log_file}: {e}")

            if folder_results:
                max_case = max(folder_results, key=lambda x: x[1])
                avg_case = max(folder_results, key=lambda x: x[2])
                logging.info("")
                logging.info("Case with highest Max value:")
                logging.info(f"{max_case[0]:<20}\tMax={max_case[1]:.1f}, Avg={max_case[2]:.1f}")
                logging.info("Case with highest Avg value:")
                logging.info(f"{avg_case[0]:<20}\tMax={avg_case[1]:.1f}, Avg={avg_case[2]:.1f}")

if __name__ == "__main__":
    main()