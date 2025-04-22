# Bellmouth 2D Simulation

This repository contains a parametric study of bellmouth inlet geometries in a 2D computational fluid dynamics (CFD) simulation using OpenFOAM. It automates the generation of geometry, meshing, simulation execution, and post-processing to investigate the aerodynamic performance of different bellmouth shapes.

---

## Acknowledgements

The WindShape team for context and application guidance

## Table of Contents

- [About the Project](#about-the-project)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## About the Project

Bellmouths are often used at the inlet of wind tunnels and ducts to smooth the airflow and reduce pressure losses. This project provides a framework to evaluate the effect of different bellmouth geometries on flow characteristics in 2D using OpenFOAM and Gmsh.

The project performs a sweep over several geometric parameters, generating and simulating multiple configurations in an automated pipeline. The goal is to help determine the influence of shape on pressure drop and velocity distribution.

Key features include:

- Parametric geometry and mesh generation using Gmsh and Python
- WSL-based automation for OpenFOAM from a Windows host
- Post-processing with integrated pressure and velocity data extraction
- Logging and reproducibility built-in for all runs

---

## Getting Started

### Prerequisites

- [Python ≥ 3.8](https://www.python.org/downloads/)
- [OpenFOAM v2412](https://openfoam.org/)
- [Gmsh](https://gmsh.info/)
- [WSL (Windows Subsystem for Linux)](https://docs.microsoft.com/en-us/windows/wsl/)
- Bash shell (via WSL)
- Optional: [Anaconda/Miniconda](https://docs.conda.io/en/latest/) for environment management

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/bellmouth-2d-simulation.git
   cd bellmouth-2d-simulation
   ```
   
2. Set up a Python environment (optional but recommended):
   
   ```bash
   conda create -n windshape python=3.10
   conda activate windshape
   ```

3. Install dependencies:
   
   ```bash
   pip install -r requirements.txt
   ```

4. Confirm OpenFOAM is installed and available in WSL:

   ```bash
   wsl
   source /usr/lib/openfoam/openfoam2412/etc/bashrc
   simpleFoam -help
   ```

## Usage
From Windows terminal (PowerShell or CMD), run:

   ```bash
   python run.py
   ```

This script will:

- Generate geometry and mesh using Gmsh
- Copy and set up OpenFOAM case directories
- Run mesh conversion, simulation, and post-processing within WSL
- Log outputs and simulation results for later analysis

You can customize geometric parameters and loop ranges in run.py.

---

## Project Structure

```php
bellmouth-2d-simulation/
├── mesh.py              # Gmsh geometry and mesh generator
├── run.py               # Main control script (simulation + post-processing)
├── ELL-case-template/        # OpenFOAM case template
├── outputs/                  # Optional directory for results
├── requirements.txt
└── README.md
```
---

## Contributing

Contributions are welcome! To propose a change:

1. Fork the repository
2. Create a new branch (```git checkout -b feature/new-feature```)
3. Commit your changes (```git commit -m 'Add new feature'```)
4. Push to the branch (```git push origin feature/new-feature```)
5. Open a Pull Request

Please ensure code adheres to project conventions and is thoroughly tested.

---

## License
This project is licensed under SOLIM's LICENSE

---

## Additional

#### Screenshots and GIFs

- Mesh and boundary conditions
- Results
- Saving data

---

#### Documentation

- [GMSH guide](https://gmsh.info/doc/texinfo/gmsh.html)
- [GMSH Python API](https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_13_1/api/gmsh.py)
