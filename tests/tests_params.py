from mesh import parse_args
import sys

def test_default_args(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['mesh.py'])
    args = parse_args()
    assert args.Mw == 12
    assert args.Kx == 0.3
