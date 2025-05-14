from utils.config import *
import gmsh

# -----------------------------
# Geometry Helper Functions
# -----------------------------

def compute_geometry_parameters(
        Mw: int, Mb: int = None, Kx: float = None, Ky: float = None
        ) -> tuple:
    """
    Compute key geometry parameters.

    Modes:
    - If only Mw is given, returns (Di,)
    - If Mb, Kx, and Ky are given, returns (Di, Db, a, b)

    Parameters:
        Mw (int): Number of WindShaper modules
        Mb (int, optional): Bellmouth reference modules
        Kx (float, optional): X-axis ellipse scale
        Ky (float, optional): Y-axis ellipse scale

    Returns:
        tuple: (Di,) or (Di, Db, a, b)
    """
    Di = Mw * MODULE_WIDTH + (Mw - 1) * MODULE_GAP  # inlet width in meters
    
    if Mb is None and Kx is None and Ky is None:
        return (Di,)
    elif Mb is not None and Kx is not None and Ky is not None:
        Db = Mb * MODULE_WIDTH + (Mb - 1) * MODULE_GAP  # bellmouth reference width in meters
        # a = 3*Kx*Db
        if Kx == 0.33:
            a_raw = Db      
        else:
            a_raw = 3 * Kx * Db
        b_raw = Ky * Db # b = Ky*Db
        
        # Ensure a >= b for ellipse definition
        a = max(a_raw, b_raw)
        b = min(a_raw, b_raw)
        return (Di, Db, a, b)
    else:
        raise ValueError("If Mb is provided, Kx and Ky must also be specified.")

 
def add_point(x: float, y:float, z:float = 0) -> int:
    """Add a point in the OCC geometry kernel."""
    return gmsh.model.occ.addPoint(x, y, z)


def add_line(p1:int, p2:int) -> int:
    """Add a straight line between two points."""
    return gmsh.model.occ.addLine(p1, p2)

def add_circle_arc(p1:int, center:int, p2:int) -> int:
    """Add a circular arc from p1 to p2 around center."""
    return gmsh.model.occ.addCircleArc(p1, center, p2)


def add_ellipse_arc(p1:int, center:int, major:int, p2:int) -> int:
    """Add an elliptical arc from p1 to p2 using major axis vector."""
    return gmsh.model.occ.addEllipseArc(p1, center, major, p2)
