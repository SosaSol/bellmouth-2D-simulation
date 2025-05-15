import pytest
from utils.geometry_utils import compute_geometry_parameters
import sys
import types
import importlib

# Mock config values (as they are imported with *)

mock_config = types.ModuleType("utils.config")
mock_config.MODULE_WIDTH = 0.5
mock_config.MODULE_GAP = 0.1
sys.modules["utils.config"] = mock_config

# Re-import after mocking config
import utils.geometry_utils as geometry_utils
importlib.reload(geometry_utils)

def test_only_Mw():
    # Di = Mw * MODULE_WIDTH + (Mw - 1) * MODULE_GAP
    Mw = 4
    expected_Di = 4 * 0.5 + (4 - 1) * 0.1  # 2.0 + 0.3 = 2.3
    result = geometry_utils.compute_geometry_parameters(Mw)
    assert result == (expected_Di,)

def test_Mw_Mb_Kx_Ky_Kx_033():
    Mw = 3
    Mb = 2
    Kx = 0.33
    Ky = 0.5
    Di = Mw * 0.5 + (Mw - 1) * 0.1  # 1.5 + 0.2 = 1.7
    Db = Mb * 0.5 + (Mb - 1) * 0.1  # 1.0 + 0.1 = 1.1
    a_raw = Db  # since Kx == 0.33
    b_raw = Ky * Db  # 0.5 * 1.1 = 0.55
    a = max(a_raw, b_raw)  # 1.1
    b = min(a_raw, b_raw)  # 0.55
    result = geometry_utils.compute_geometry_parameters(Mw, Mb, Kx, Ky)
    assert result == (Di, Db, a, b)

def test_Mw_Mb_Kx_Ky_Kx_not_033():
    Mw = 2
    Mb = 3
    Kx = 0.4
    Ky = 0.6
    Di = Mw * 0.5 + (Mw - 1) * 0.1  # 1.0 + 0.1 = 1.1
    Db = Mb * 0.5 + (Mb - 1) * 0.1  # 1.5 + 0.2 = 1.7
    a_raw = 3 * Kx * Db  # 3 * 0.4 * 1.7 = 2.04
    b_raw = Ky * Db      # 0.6 * 1.7 = 1.02
    a = max(a_raw, b_raw)  # 2.04
    b = min(a_raw, b_raw)  # 1.02
    result = geometry_utils.compute_geometry_parameters(Mw, Mb, Kx, Ky)
    assert result == (Di, Db, a, b)

def test_missing_Kx_raises():
    with pytest.raises(ValueError):
        geometry_utils.compute_geometry_parameters(2, Mb=2, Kx=None, Ky=0.5)

def test_missing_Ky_raises():
    with pytest.raises(ValueError):
        geometry_utils.compute_geometry_parameters(2, Mb=2, Kx=0.4, Ky=None)

def test_missing_Mb_raises():
    with pytest.raises(ValueError):
        geometry_utils.compute_geometry_parameters(2, Mb=None, Kx=0.4, Ky=0.5)

def test_missing_Mw_raises():
    with pytest.raises(TypeError):
        geometry_utils.compute_geometry_parameters(Mw=None, Mb=2, Kx=0.4, Ky=0.5)

def test_missing_Mw_and_Mb_raises():
    with pytest.raises(TypeError):
        geometry_utils.compute_geometry_parameters(Mw=None, Mb=None, Kx=0.4, Ky=0.5)    
