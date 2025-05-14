"""
Configuration file for bellmouth-2D-simulation.

Defines project constants and parameters for geometry and logging.
"""

import logging
import math
from pathlib import Path

# =========================
# Project Geometry Constants
# =========================

# Module dimensions (in meters)
MODULE_WIDTH = 240.0e-3  # Width of the module

# Gap dimensions (in meters)
MODULE_GAP = 2.5e-3             # Primary gap
GAP = 42.56e-3          # Secondary gap

# Lengths (in meters)
L1 = 168.66e-3           # First section length
L2 = 130.29e-3           # Second section length
L3 = 39.86e-3            # Third section length
L = L1 + L2 + L3         # Total length

# Heights (in meters)
H1 = 199.94e-3           # First section height
# H1 = 217.5e-3          # Alternative height (commented out)
H2 = 238.94e-3           # Second section height

# Angles (in degrees)
ANGLE = 180 - 171.79     # Angle for geometry calculation

# Derived values
DELTA_H = L2 * math.tan(math.radians(ANGLE))  # Height difference due to angle

# =========================
# Width Dimensions (in meters)
# =========================

t_22 = 55e-3     # For 2x2 to 3x3 configurations
t_44 = 75e-3     # For 4x4 to 12x12 configurations
t_tilt = 125e-3  # For tilting machine
t_rack = 528e-3  # With racks info and PDBox

# =========================
# Ellipse Parameter (in meters)
# =========================

a_ell = 70e-3    # Ellipse semi-major axis

# =========================
# Output Directory
# =========================

SAVE_ROOT = Path.cwd()   # Current working directory as save root

# =========================
# Logging Configuration
# =========================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s: %(message)s"
)

# =========================
# End of config.py
# =========================