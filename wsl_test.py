import subprocess
from pathlib import Path
import time

# --- Configuration ---
PYTHON_EXECUTABLE = Path(r"C:/Users/solim/miniconda3/envs/windshape/python.exe")
SCRIPT_PATH = Path("mesh_test.py")

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
WSL_BASE = "/home/solim/OpenFOAM/solim-v2412/run/bellmouth_2D"
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
        print(f"[ERROR] Command failed: {command}")
        print(result.stderr)
    else:
        print(result.stdout)

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
        str(PYTHON_EXECUTABLE),
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
        print(f"[ERROR] Mesh generation failed for {fname}")
        print(f"[DETAILS] {e}")
    else:
        print(f"[INFO] Mesh generation successful for {fname}")

def main():
    total_jobs = len(MW_RANGE) * len(MB_RANGE)
    job_counter = 0
    zero_time = time.time()

    for Mw in MW_RANGE:
        for Mb in MB_RANGE:
            job_counter += 1
            fname = f"ELL-{Mw}-{Mb}-{int(Kx*100)}-{int(Ky*100)}-{int(r*1e3)}-{int(t*1e3)}"
            case_path = wsl_path(fname)

            print(f"\n{'-'*60}")
            print(f"JOB {job_counter} of {total_jobs}: {fname}")
            print(f"{'-'*60}")

            # 1. Copy template case if it doesn't exist
            if run_command_wsl(f"test -d {case_path}"):
                print(f"[SKIP] Case folder already exists: {case_path}")
            else:
                print(f"[INFO] Copying template case to: {case_path}")
                run_command_wsl(f"cp -r {WSL_BASE}/{CASE_TEMPLATE} {case_path}")
                run_command_wsl(f"chmod +x Allrun", cwd=case_path)

            # 2. Check mesh
            if mesh_exists(fname):
                print(f"[SKIP] Mesh already exists for {fname}, skipping mesh generation.")
            else:
                print(f"[INFO] Generating mesh for {fname}")
                run_gmsh_script(Mw, Mb, fname)

            # 3. Run Allrun to execute the full simulation
            print(f"[INFO] Running OpenFOAM simulation via Allrun script")
            run_command_wsl("./Allrun", cwd=case_path)

            print(f"[DONE] Job {job_counter} of {total_jobs} completed: {fname}")
            print(f"{'-'*60}")

    total_time = time.time() - zero_time
    print(f"\n All {total_jobs} jobs completed in {total_time:.2f} seconds")

if __name__ == "__main__":
    main()
