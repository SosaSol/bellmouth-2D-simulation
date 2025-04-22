import pytest
import sys
import os

# Add the path to your mesh.py file
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from mesh import compute_geometry_parameters

def test_compute_geometry_parameters_basic():
    Di, Db, a, b = compute_geometry_parameters(Mw=10, Mb=15, Kx=0.3, Ky=0.5)
    assert Di > 0
    assert Db > 0
    assert a >= b
    assert isinstance(Di, float)
