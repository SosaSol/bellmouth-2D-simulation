import subprocess
from pathlib import Path
import time
import logging

# --- Configuration ---
# PYTHON_EXECUTABLE = Path(r"C:/Users/solim/miniconda3/envs/windshape/python.exe")
SCRIPT_PATH = Path("mesh.py")

# --- logging configuration ---
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s")

# Simulation constants
Kx, Ky = 0.3, 0.3
r, t = 10e-3, 5e-3
L = 0.3
xmin, ymin = 0, 0
xmax, ymax = 25, 25
nt = 15

# Loop range
MW_RANGE = range(2, 13)
MB_RANGE = range(2, 13)

# WSL Paths
WSL_BASE = str(Path(__file__).parent.resolve())
CASE_TEMPLATE = "ELL-case-template"

def wsl_path(*parts):
    """Join and return WSL path"""
    return f"{WSL_BASE}/{'/'.join(parts)}"

import subprocess

def run_command_wsl(command, cwd=None):
    """Run a command in WSL with OpenFOAM environment sourced."""
    source_cmd = "source /usr/lib/openfoam/openfoam2412/etc/bashrc"
    cd_cmd = f"cd {cwd} && " if cwd else ""
    full_cmd = f"{source_cmd} && {cd_cmd}{command}"

    result = subprocess.run(
        ["wsl", "bash", "-c", full_cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # Determine if it was successful
    if result.returncode != 0:
        # Optionally suppress expected 'failures' for checks
        if 'test -d' in command or 'ls' in command:
            return False
        logging.info(f"Command failed: {command}")
        logging.info(result.stderr)
    else:
        logging.info(result.stdout)

    return result


def mesh_exists(fname):
    """Check if the mesh file already exists in WSL."""
    trisurf_dir = wsl_path(fname, "constant", "triSurface")
    check_command = f"ls {trisurf_dir}/*.msh 2>/dev/null"
    result = run_command_wsl(check_command)
    return result

def run_gmsh_script(Mw, Mb, fname):
    """Run the mesh generation script with specified parameters."""
    command = [
        str("python"),
        str(SCRIPT_PATH),
        "--Mw", str(Mw),
        "--Mb", str(Mb),
        "--Kx", str(Kx),
        "--Ky", str(Ky),
        "--r", str(r),
        "--t", str(t),
        "--L", str(L),
        "--xmin", str(xmin),
        "--ymin", str(ymin),
        "--xmax", str(xmax),
        "--ymax", str(ymax),
        "--nt", str(nt),
    ]

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logging.info(f"Mesh generation failed for {fname}")
        logging.info(f"{e}")
    else:
        logging.info(f"Mesh generation successful for {fname}")

def main():
    total_jobs = len(MW_RANGE) * len(MB_RANGE)
    job_counter = 0
    zero_time = time.time()

    for Mw in MW_RANGE:
        for Mb in MB_RANGE:
            job_counter += 1
            fname = f"ELL-{Mw}-{Mb}-{int(Kx*100)}-{int(Ky*100)}-{int(r*1e3)}-{int(t*1e3)}"
            case_path = wsl_path(fname)

            logging.info(f"\n{'-'*60}")
            logging.info(f"JOB {job_counter} of {total_jobs}: {fname}")
            logging.info(f"{'-'*60}")

            # 1. Copy template case if it doesn't exist
            if run_command_wsl(f"test -d {case_path}"):
                logging.info(f"Case folder already exists: {case_path}")
            else:
                logging.info(f"Copying template case to: {case_path}")
                run_command_wsl(f"cp -r {WSL_BASE}/{CASE_TEMPLATE} {case_path}")
                run_command_wsl(f"chmod +x Allrun", cwd=case_path)

            # 2. Check mesh
            if mesh_exists(fname):
                logging.info(f"Mesh already exists for {fname}, skipping mesh generation.")
            else:
                logging.info(f"Generating mesh for {fname}")
                run_gmsh_script(Mw, Mb, fname)

            # 3. Run Allrun to execute the full simulation
            logging.info(f"Running OpenFOAM simulation via Allrun script")
            run_command_wsl("./Allrun", cwd=case_path)

            logging.info(f"[DONE] Job {job_counter} of {total_jobs} completed: {fname}")
            logging.info(f"{'-'*60}")

    total_time = time.time() - zero_time
    logging.info(f"\n All {total_jobs} jobs completed in {total_time:.2f} seconds")

if __name__ == "__main__":
    main()
