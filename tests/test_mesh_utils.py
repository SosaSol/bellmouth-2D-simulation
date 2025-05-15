import math
import pytest
from utils.mesh_utils import boundaryLayerParameters

def test_boundary_layer_parameters_typical():
    # Typical values for air at room temperature and moderate velocity
    U_inf = 16.0
    nu = 15.06e-6
    x = 0.3
    y_plus = 0.95
    n_layers = 20
    yp_factor = 1.0

    y1, delta99, expansion_ratio = boundaryLayerParameters(
        U_inf=U_inf, nu=nu, x=x, y_plus=y_plus, n_layers=n_layers, yp_factor=yp_factor
    )

    assert y1 > 0
    assert delta99 > 0
    assert 1.0 < expansion_ratio < 2.5  # Reasonable geometric growth ratio

def test_boundary_layer_parameters_low_layers_warns(caplog):
    # n_layers < 5 should trigger a warning
    with caplog.at_level("WARNING"):
        y1, delta99, expansion_ratio = boundaryLayerParameters(
            U_inf=10.0, nu=1.5e-5, x=0.2, y_plus=1.0, n_layers=3
        )
        assert "n_layers < 5" in caplog.text
        assert y1 > 0
        assert delta99 > 0

def test_boundary_layer_parameters_yp_factor_scaling():
    # yp_factor should inversely scale y1, but not affect delta99 or expansion_ratio significantly
    params_default = boundaryLayerParameters(
        U_inf=10, nu=1.5e-5, x=0.2, y_plus=1.0, n_layers=10, yp_factor=1.0
    )
    params_half = boundaryLayerParameters(
        U_inf=10, nu=1.5e-5, x=0.2, y_plus=1.0, n_layers=10, yp_factor=0.5
    )
    y1_default, delta99_default, ratio_default = params_default
    y1_half, delta99_half, ratio_half = params_half

    # y1 should half when yp_factor is halved
    assert math.isclose(y1_half*2, y1_default, rel_tol=1e-6)
    # delta99 should remain unchanged
    assert math.isclose(delta99_half, delta99_default, rel_tol=1e-6)

def test_boundary_layer_parameters_high_reynolds():
    # For high Re_x, delta99 formula changes
    U_inf = 100.0
    nu = 1.5e-5
    x = 2.0
    y_plus = 1.0
    n_layers = 15
    y1, delta99, expansion_ratio = boundaryLayerParameters(
        U_inf=U_inf, nu=nu, x=x, y_plus=y_plus, n_layers=n_layers
    )
    assert y1 > 0
    assert delta99 > 0
    assert 1.0 < expansion_ratio < 2.5
