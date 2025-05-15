# Bellmouth 2D Simulation

This repository provides a fully automated workflow for parametric CFD studies of 2D bellmouth inlet geometries using OpenFOAM and Gmsh. The project enables users to generate, mesh, simulate, and post-process multiple bellmouth configurations efficiently, supporting rapid design exploration and analysis.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [References](#references)

---

## Overview

Bellmouth inlets are critical for minimizing pressure losses and ensuring smooth airflow in wind tunnels and duct systems. This project automates the process of evaluating how different 2D bellmouth shapes affect flow characteristics such as pressure drop and velocity distribution.

The workflow includes:

- Parametric geometry and mesh generation with Gmsh (via Python)
- Automated OpenFOAM case setup and execution in WSL
- Post-processing for extracting and logging key aerodynamic metrics

---

## Features

- **Parametric Study:** Easily configure and sweep geometric parameters for bellmouth shapes.
- **Automated Pipeline:** One-command execution for geometry, meshing, simulation, and post-processing.
- **Cross-Platform:** Designed for Windows hosts using WSL for Linux-based OpenFOAM.
- **Reproducibility:** All runs are logged with configuration and results for traceability.
- **Extensible:** Modular Python scripts for geometry, simulation control, and data extraction.

---

## Requirements

- Windows 10/11 with [WSL](https://docs.microsoft.com/en-us/windows/wsl/) (Ubuntu recommended)
- [Python ≥ 3.8](https://www.python.org/downloads/)
- [OpenFOAM v2412](https://openfoam.org/) (installed in WSL)
- [Gmsh](https://gmsh.info/) (installed in WSL)
- Bash shell (via WSL)
- Optional: [Anaconda/Miniconda](https://docs.conda.io/en/latest/) for Python environment management
- Additional Linux packages: `libglu1-mesa`, `libxft2`

---

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/SosaSol/bellmouth-2d-simulation.git
   cd bellmouth-2d-simulation
   ```

2. **Set up Python environment (optional but recommended):**

   ```bash
   conda create -n bellmouth python=3.10
   conda activate bellmouth
   ```

3. **Install Python dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

4. **Install OpenFOAM and Gmsh in WSL:**
   - Follow [OpenFOAM installation guide](https://openfoam.org/download/)
   - Install Gmsh: `sudo apt install gmsh`
   - Install required libraries: `sudo apt install libglu1-mesa libxft2`

---

## Usage

1. **Configure parameters:**  
   Edit `utils/config.py` to set geometric parameter ranges and simulation options.  
   Additionally, use `run_all.sh` to configure specific geometry settings and define parameter sweep options before running the workflow.

2. **Run the workflow from wsl terminal:**

   ```bash
   chmod +x run_all.sh
   ./run_all.sh
   ```

   This will:

   - Generate geometry and mesh for each configuration
   - Set up and run OpenFOAM simulations
   - Post-process results and save logs/outputs

3. **Results:**  
   Output data and logs are saved in the `outputs/` directory for analysis.

---

## Project Structure

```text
bellmouth-2d-simulation/
├── mesh.py                          # Meshing script for parametric and curved bellmouth geometries
├── mesh_straight.py                 # Meshing script for straight inlet geometries
├── mesh_advanced_no_bellmouth.py    # Meshing script for cases without bellmouth features
├── mesh_advanced_with_bellmouth.py  # Advanced meshing for curved bellmouth inlets
├── mesh_aeropack.py                 # Meshing for aerodynamic package configurations
├── utils/
│   ├── config.py                    # Parameter configuration for geometry and simulation
│   ├── geometry.py                  # Parametric geometry generation functions
│   ├── meshing.py                   # Mesh generation utilities using Gmsh
│   └── ...                          # Additional utility scripts
├── tests/
│   └── ...
├── ELL-case-template/               # OpenFOAM case template directory
├── outputs/                         # Simulation results and logs
├── postProcessingOutputs/           # Post-processed results (plots, processed data)
├── plotResiduals.py                 # Plotting simulation residuals
├── plotLiveResiduals.py             # Compare results across simulation cases
├── run_all.sh                       # Automates the full workflow
├── clean_all.sh                     # Cleans project folders and outputs
├── requirements.txt                 # Python dependencies
├── updateBoundaryConditions.py      # Modify OpenFOAM boundary conditions
└── README.md
```

---

## Contributing

Contributions are welcome!  
To propose changes:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit and push your changes
4. Open a Pull Request

Please follow project conventions and test your code.

---

## License

This project is licensed under SOLIM's LICENSE. See `LICENSE` file for details.

---

## Acknowledgements

- Tony G.
- Aurélien W.
- Sergio M.
- The whole WindShape team for project context and application guidance

---

## References

- [Gmsh Documentation](https://gmsh.info/doc/texinfo/gmsh.html)
- [Gmsh Python API](https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/api/gmsh.py)
- [OpenFOAM User Guide](https://cfd.direct/openfoam/user-guide/)
- [WSL Documentation](https://docs.microsoft.com/en-us/windows/wsl/)
