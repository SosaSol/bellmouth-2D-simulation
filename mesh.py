# standard library
import os
import sys
import argparse
import multiprocessing
from pathlib import Path
import logging
import math
from scipy.optimize import fsolve

# third‑party
import gmsh

# project constants
MM = 1e-3
MODULE_WIDTH = 240.0e-3  # m
GAP = 2.5e-3  # m
OPENFOAM_BASHRC = "/usr/lib/openfoam/openfoam2412/etc/bashrc"

SAVE_ROOT = Path.cwd()

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s")

# -----------------------------
# Parser
# -----------------------------

def parse_args():
    """
    Parse command line arguments for mesh generation.
    """
    parser = argparse.ArgumentParser(description="Generate GMSH mesh for inlet geometry.")
    
    parser.add_argument("--Mw", type=int, default=12, help="Number of wires")
    parser.add_argument("--Mb", type=int, default=12, help="Bellmouth multiplier")
    
    parser.add_argument("--Kx", type=float, default=0.3, help="X-axis ellipse scale")
    parser.add_argument("--Ky", type=float, default=0.3, help="Y-axis ellipse scale")
    
    parser.add_argument("--r", type=float, default=10e-3, help="Radius of curvature (in meters)")
    parser.add_argument("--t", type=float, default=5e-3, help="Wall thickness (in meters)")

    parser.add_argument("--L", type=float, default=0.3, help="Length of straight section (in meters)")

    parser.add_argument("--xmin", type=float, default=0, help="X min (in meters)")
    parser.add_argument("--ymin", type=float, default=0, help="Y min (in meters)")
    parser.add_argument("--xmax", type=float, default=25, help="X max (in meters)")
    parser.add_argument("--ymax", type=float, default=25, help="Y max (in meters)")

    parser.add_argument("--nt", type=int, default=10, help="Number of threads for GMSH")
    parser.add_argument("--sd", type=str, default=None, help="Save directory")
    return parser.parse_args()

# -----------------------------
# Geometry Helper Functions
# -----------------------------
def compute_geometry_parameters(
        Mw: int, Mb: int, Kx: float, Ky: float
        ) -> tuple[float, float, float, float]:
    """
    Compute key geometry parameters based on input multipliers.
    Returns:
        Di: inlet width
        Db: bellmouth reference width
        a: semi-major axis (half bellmouth width)
        b: semi-minor axis (half bellmouth height)
    """
    # Di and Db are calculated as:
    Di = Mw * MODULE_WIDTH + (Mw - 1) * GAP  # inlet width in meters
    Db = Mb * MODULE_WIDTH + (Mb - 1) * GAP  # bellmouth reference width in meters

    # Ensure a >= b for ellipse definition
    a_raw = Db / 2
    b_raw = Ky * Db / 2
    a = max(a_raw, b_raw)
    b = min(a_raw, b_raw)
    return Di, Db, a, b


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

# -----------------------------
# Build 2D Geometry
# -----------------------------
def create_geometry(
        Mw:int, Mb:int, Kx: float, Ky:float, 
        r:float, t:float, L:float, 
        xmin:float, ymin:float, xmax:float, ymax:float
        ) -> tuple[int, list[int], list[int]]:
    """
    Build the parametric 2D bellmouth inlet geometry.
    Returns:
        surface: the tag of the created surface
        curve_tags: list of all curve tags in order
        edge_points: two point tags for boundary layer 'fan'
    """
    logging.info("Building geometry...")

    # Compute parameters
    Di, Db, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)
    y_inlet = ymin + Di / 2

    # Key x positions
    x1 = xmax - L
    x2 = x1 - a

    # Centers for arcs
    ell_center    = add_point(x1, y_inlet + b)
    circ_center   = add_point(x2 + r, y_inlet + b)
    circ_center2  = add_point(x2 + r, y_inlet + b + r - t/2)

    # Create boundary points
    pts = []
    pts.append(add_point(xmin, ymin))            # p1
    pts.append(add_point(xmax, ymin))            # p2
    pts.append(add_point(xmax, y_inlet))         # p3
    pts.append(add_point(x1, y_inlet))           # p4
    pts.append(add_point(x2, y_inlet + b))       # p5
    pts.append(add_point(x2 + r, y_inlet + b + r))           # p6
    pts.append(add_point(x2 + r, y_inlet + b + r - t))       # p7
    pts.append(add_point(x2 + t, y_inlet + b))               # p8
    pts.append(add_point(x1, y_inlet + t))       # p9
    pts.append(add_point(xmax, y_inlet + t))     # p10
    pts.append(add_point(xmax, ymax))            # p11
    pts.append(add_point(xmin, ymax))            # p12

    # Store edge points for boundary layer
    edge_points = [pts[2], pts[9]]  # p3, p10

    # Build curves sequentially
    curve_tags = []
    curve_tags.append(add_line(pts[0], pts[1]))           # bottom
    curve_tags.append(add_line(pts[1], pts[2]))           # inlet lower wall
    curve_tags.append(add_line(pts[2], pts[3]))           # straight before ellipse
    curve_tags.append(add_ellipse_arc(pts[3], ell_center, pts[4], pts[4]))
    curve_tags.append(add_circle_arc (pts[4], circ_center,  pts[5]))
    curve_tags.append(add_circle_arc (pts[5], circ_center2, pts[6]))
    curve_tags.append(add_circle_arc (pts[6], circ_center,  pts[7]))
    curve_tags.append(add_ellipse_arc(pts[7], ell_center, pts[7], pts[8]))
    curve_tags.append(add_line(pts[8], pts[9]))           # inlet upper wall
    curve_tags.append(add_line(pts[9], pts[10]))          # right vertical
    curve_tags.append(add_line(pts[10], pts[11]))         # top horizontal
    curve_tags.append(add_line(pts[11], pts[0]))          # left vertical

    # Create surface
    loop_tag = gmsh.model.occ.addCurveLoop(curve_tags)
    surface = gmsh.model.occ.addPlaneSurface([loop_tag])
    gmsh.model.occ.synchronize()

    return surface, curve_tags, edge_points

# -----------------------------
# Boundary Layer
# -----------------------------

def layers_needed_numeric(y1, delta, r):
    def func(N):
        return y1 * (1 - r**N) / (1 - r) - delta
    N_guess = 20
    N_sol = fsolve(func, N_guess)
    return math.ceil(N_sol[0])

def compute_boundary_layer_info(U_inf: float, nu: float, x:float, y_plus:float=0.95, expansion_ratio:float=1.05):
    """
    Computes boundary layer parameters for a turbulent boundary layer over a flat plate.
    
    Parameters:
    - U_inf: Freestream velocity [m/s]
    - nu: Kinematic viscosity [m^2/s]
    - x: Distance from leading edge [m] (for estimating Cf and delta99)
    - y_plus: Desired dimensionless first cell height
    - growth_rate: Geometric growth rate between layers
    - n_layers: Number of boundary layer layers
    
    Returns:
    - Dictionary with y1, delta99, total height, and layer details
    """
        
    if expansion_ratio <= 1.0:
        raise ValueError("Expansion ratio must be > 1.0")

    
    Re_x = U_inf * x / nu               # Reynolds number based on x 
    Cf = 0.026 / Re_x**0.2              # Skin friction coefficient for turbulent boundary layer (empirical)
    u_tau = U_inf * math.sqrt(Cf / 2)   # Friction velocity
    y1 = y_plus * nu / u_tau            # First cell height for desired y+
    delta99 = 0.37 * x / Re_x**0.2      # Boundary layer thickness
    
    # Compute required number of layers
    n_layers = math.log(1 + ((expansion_ratio - 1) * delta99 / y1)) / math.log(expansion_ratio)
    n_layers = math.ceil(n_layers)          # Round up to nearest integer

    return {
        "Re_x": Re_x,
        "Cf": Cf,
        "u_tau": u_tau,
        "y1": y1,
        "delta99": delta99,
        "n_layers": n_layers,
        "expansion_ratio": expansion_ratio
    }

def apply_boundary_layer(curve_tags:list, edge_points:list, bl_thickness:float=3e-3):
    """
    Apply a structured boundary layer on selected curves and fan points.
    """
    bl = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(bl, "CurvesList", curve_tags)
    gmsh.model.mesh.field.setNumbers(bl, "PointsList", edge_points)
    gmsh.model.mesh.field.setNumber(bl, "Thickness", bl_thickness)
    gmsh.model.mesh.field.setNumber(bl, "Size", 5e-5)
    gmsh.model.mesh.field.setNumber(bl, "Ratio", 1.05)
    gmsh.model.mesh.field.setNumber(bl, "NbLayers", 30)
    gmsh.model.mesh.field.setNumber(bl, "Quads", 1)
    gmsh.model.mesh.field.setAsBoundaryLayer(bl)
    gmsh.model.occ.synchronize()

# -----------------------------
# Mesh Refinement
# -----------------------------
def refine_mesh(curve_tags:list, xmin:float, xmax:float, ymin:float, ymax:float):
    """
    Add distance + threshold + box fields and combine for smooth grading.
    """
    # Distance from bellmouth wall
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "EdgesList", curve_tags[2:9])
    gmsh.model.mesh.field.setNumber(dist, "Sampling", 100)

    # Threshold based on distance
    thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thr, "InField", dist)
    gmsh.model.mesh.field.setNumber(thr, "SizeMin", 1e-3)
    gmsh.model.mesh.field.setNumber(thr, "SizeMax", 2)
    gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.02)
    gmsh.model.mesh.field.setNumber(thr, "DistMax", 0.03)

    # Box refinement near inlet corner
    box = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(box, "VIn", 1e-2)
    gmsh.model.mesh.field.setNumber(box, "VOut", 1)
    gmsh.model.mesh.field.setNumber(box, "XMin", xmax-3)
    gmsh.model.mesh.field.setNumber(box, "XMax", xmax)
    gmsh.model.mesh.field.setNumber(box, "YMin", ymin)
    gmsh.model.mesh.field.setNumber(box, "YMax", ymin+3)
    gmsh.model.mesh.field.setNumber(box, "ZMin", -2)
    gmsh.model.mesh.field.setNumber(box, "ZMax",  2)
    gmsh.model.mesh.field.setNumber(box, "Thickness", 10)

    # Combine fields
    mfield = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mfield, "FieldsList", [thr, box])
    gmsh.model.mesh.field.setAsBackgroundMesh(mfield)
    gmsh.model.occ.synchronize()

# -----------------------------
# Extrusion and Physical Groups
# -----------------------------
def extrude_and_group(surface:int):
    """
    Extrude the 2D surface into a thin 3D volume by one element,
    then define physical groups for CFD boundaries and volume.
    """
    # Extrude surface into volume
    extruded = gmsh.model.occ.extrude([(2, surface)], 0, 0, 1.0,
                                       numElements=[1], recombine=True)
    gmsh.model.occ.synchronize()

    inletTag        = [11, 12, 13]
    outletTag       = [3]
    frontBackTag    = [1, 14]
    wallTag         = [4, 5, 6, 7, 8, 9, 10]
    symmetryTag     = [2]
    volumeTag       = [1]

    # Create surface boundary groups
    gmsh.model.addPhysicalGroup(2, inletTag,        name="inlet")
    gmsh.model.addPhysicalGroup(2, outletTag,       name="outlet")
    gmsh.model.addPhysicalGroup(2, frontBackTag,    name="frontAndBack")
    gmsh.model.addPhysicalGroup(2, wallTag,         name="wall")
    gmsh.model.addPhysicalGroup(2, symmetryTag,     name="symmetry")

    # Create internal volume group
    gmsh.model.addPhysicalGroup(3, volumeTag,       name="internal")

# -----------------------------
# Save Mesh
# -----------------------------
def save_mesh(fname:str , save_path: Path):
    save_path.mkdir(parents=True, exist_ok=True)
    msh = save_path / f"{fname}.msh"
    gmsh.write(str(msh))
    logging.info(f"mesh saved to {msh}")

# -----------------------------
# Print Information
# -----------------------------    
def print_info(
        Mw:int, Mb:int, Kx:float, Ky:float, 
        r:float, t:float, L:float, 
        xmin:float, ymin:float, xmax:float, ymax:float, 
        nt: int, fname:str):
    """
    Print the parameters used for mesh generation.
    """
    logging.info("Geometry parameters:")
    logging.info(f"          - Mesh tag suffix            : {fname}")
    logging.info(f"          - WindShaper modules         : {Mw}")
    logging.info(f"          - Bellmouth reference modules: {Mb}")
    logging.info(f"          - Ellipse scaling (Kx, Ky)   : {Kx:.2f}, {Ky:.2f}")
    logging.info(f"          - Radius of curvature (mm)   : {r * 1e3:.1f}")
    logging.info(f"          - Wall thickness (mm)        : {t * 1e3:.1f}")
    logging.info(f"          - Straight section len (mm)  : {L * 1e3:.1f}")

    logging.info(f"Domain bounds -> X: [{xmin}, {xmax}], Y: [{ymin}, {ymax}]")
    logging.info(f"Using {nt} threads for GMSH.")


# -----------------------------
# Validate arguments
# -----------------------------
def validate_args(args):
    """
    Validate input arguments and raise ValueError for invalid entries.
    """
    if args.Mw <= 0:
        raise ValueError("Mw (Number of wires) must be greater than 0.")
    if args.Mb <= 0:
        raise ValueError("Mb (Bellmouth multiplier) must be greater than 0.")
    if not (0 < args.Kx <= 1):
        raise ValueError("Kx must be in the range (0, 1].")
    if not (0 < args.Ky <= 1):
        raise ValueError("Ky must be in the range (0, 1].")
    if args.r <= 0:
        raise ValueError("Radius r must be positive.")
    if args.t <= 0:
        raise ValueError("Wall thickness t must be positive.")
    if args.L <= 0:
        raise ValueError("Length L must be positive.")
    if args.xmin >= args.xmax:
        raise ValueError("xmin must be less than xmax.")
    if args.ymin >= args.ymax:
        raise ValueError("ymin must be less than ymax.")
    if args.nt <= 0:
        raise ValueError("Number of threads (nt) must be positive.")


# -----------------------------
# Main Execution
# -----------------------------
def main(Mw:int=12, Mb:int=12, Kx:float=0.33, Ky:float=0.33, 
         r:float=10e-3, t:float=5e-3, L:float=0.3, 
         xmin:float=0, ymin:float=0, xmax:float=25, ymax:float=25,
        nt:int=10, sd:str=None):
    """
    Main function to generate the mesh using GMSH.
    """

    # Define name and path for mesh file
    fname = f"ELL-{Mw}-{Mb}-{int(Kx*100)}-{int(Ky*100)}-{int(r*1e3)}-{int(t*1e3)}"
    save_path = Path(sd) if sd else Path("meshes")
    # logging.info(f"Save path: {save_path}")
    

    # Initialize and add model
    gmsh.initialize()
    gmsh.model.add(fname)
    # num_threads = min(multiprocessing.cpu_count() // 2, nt)
    gmsh.option.setNumber("General.NumThreads", nt)     # Set number of threads for GMSH
    
    # Print info
    print_info(Mw, Mb, Kx, Ky, r, t, L, xmin, ymin, xmax, ymax, nt, fname)
    # Create geometry
    surf, curves, epts = create_geometry(Mw, Mb, Kx, Ky, r, t, L, xmin, ymin, xmax, ymax)
    # Apply boundary layer
    apply_boundary_layer(curve_tags=curves[2:9], edge_points=epts)
    # Refine mesh
    refine_mesh(curve_tags=curves, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

    # Generate 2D mesh
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.model.mesh.generate(2)

    # Extrude and define physical groups
    extrude_and_group(surface=surf)
    
    # Generate 3D mesh and save
    gmsh.model.mesh.generate(3)
    save_mesh(fname, save_path=save_path)

    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()

    gmsh.finalize()

# -----------------------------
# Entry point for script execution
def run():
    args = parse_args()
    try:
        validate_args(args)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(1)
    
    main(args.Mw, args.Mb, args.Kx, args.Ky, args.r, args.t,
         args.L, args.xmin, args.ymin, args.xmax, args.ymax,
         args.nt, args.sd)

if __name__ == "__main__":
    run()