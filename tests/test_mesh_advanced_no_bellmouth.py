import math
import pytest
from mesh_advanced_no_bellmouth import compute_geometry_parameters
from .utils.config import MODULE_WIDTH, MODULE_GAP

@pytest.mark.parametrize("Mw,expected", [
    (1, MODULE_WIDTH),
    (2, 2 * MODULE_WIDTH + MODULE_GAP),
    (3, 3 * MODULE_WIDTH + 2 * MODULE_GAP),
    (10, 10 * MODULE_WIDTH + 9 * MODULE_GAP),
])
def test_compute_geometry_parameters(Mw, expected):
    result = compute_geometry_parameters(Mw)
    assert math.isclose(result, expected, rel_tol=1e-9)

def test_compute_geometry_parameters_zero():
    with pytest.raises(TypeError):
        compute_geometry_parameters(0)

def test_compute_geometry_parameters_negative():
    with pytest.raises(TypeError):
        compute_geometry_parameters(-5)

