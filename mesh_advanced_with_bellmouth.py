# standard library
import sys
import argparse
from pathlib import Path
import logging
import math
from scipy.optimize import newton

# third‑party
import gmsh

# project constants
MODULE_WIDTH = 240.0e-3  # m
GAP = 2.5e-3  # m

L1 = 168.66e-3  # m
L2 = 130.29e-3  # m
L3 = 39.86e-3   # m
L = L1 + L2 + L3 # m

H1 = 199.94e-3  # m
# H1 = 217.5e-3  # m
H2 = 238.94e-3  # m

ANGLE = 180 - 171.79  # degrees
GAP2 = 42.56e-3  # m
DELTA_H = L2 * math.sin(ANGLE * math.pi/180)  # m

t = GAP2/2

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
    
    parser.add_argument("--Kx", type=float, default=0.33, help="X-axis ellipse scale")
    parser.add_argument("--Ky", type=float, default=0.33, help="Y-axis ellipse scale")
    
    parser.add_argument("--r", type=float, default=50e-3, help="Radius of curvature (in meters)")

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
        Di: inlet width [m]
        Db: bellmouth reference width [m]
        a: semi-major axis (half bellmouth width) [m]
        b: semi-minor axis (half bellmouth height) [m]
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
        r:float, 
        xmin:float, ymin:float, xmax:float, ymax:float
        ) -> tuple[int, list[int], list[int]]:
    """
    Build the parametric 2D bellmouth inlet geometry.
    Returns:
        surface: the tag of the created surface
        curve_tags: list of all curve tags in order
        edge_points: point tags for boundary layer end
    """
    logging.info("Building geometry...")

    # Parametric Inlet Geometry

    # Key x positions

    x1 = xmax - L
    x2 = xmax - L2 - L3
    x3 = xmax - L3

    # Create boundary points
    pts = []
    pts.append(add_point(xmin, ymin))            # p1

    # Store edge points for boundary layer
    edge_points = []

    # Store corner points for refinement
    corner_points = []

    # If Mw odd
    if Mw % 2 == 1:
        pts.append(add_point(xmax, ymin))            # p2         
        N = (Mw + 1) // 2
        for i in range(N):
            y1 = H1 * (i+1/2) + GAP2 * i
            y2 = y1 + DELTA_H
            y4 = y1 + GAP2
            y3 = y4 - DELTA_H

            pts.append(add_point(xmax, y2))
            edge_points.append(pts[-1])  # Store edge point for boundary layer
            pts.append(add_point(x3, y2))
            pts.append(add_point(x2, y1))
            pts.append(add_point(x1, y1))
            corner_points.append(pts[-1])  # Store corner point for refinement

            if i == N-1:
                y4 = y1 + t
                y3 = y4
                Di, Db, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)

                x0 = x1 - a
                y0 = y1 + b

                # ellipse points
                pts.append(add_point(x0, y0))
                pts.append(add_point(x0 + r, y0 + r))
                pts.append(add_point(x0 + r, y0 + r -t))
                pts.append(add_point(x0 + t, y0))


            pts.append(add_point(x1, y4))
            corner_points.append(pts[-1])  # Store corner point for refinement

            pts.append(add_point(x2, y4))
            pts.append(add_point(x3, y3))
            pts.append(add_point(xmax, y3))
            edge_points.append(pts[-1])  # Store edge point for boundary layer
            
    # If Mw even
    elif Mw % 2 == 0:
        N = Mw // 2 + 1
        for i in range(N):
            y1 = H1 * i + GAP2 * i/2
            y2 = y1 + DELTA_H
            y4 = y1 + GAP2
            y3 = y4 - DELTA_H

            if i == 0:
                pts.append(add_point(x1, ymin))            # p2
                edge_points.append(pts[-1])  # Store edge point for boundary layer
                y4 = GAP2/2
                y3 = y4 - DELTA_H
            elif i!= 0:
                pts.append(add_point(xmax, y2))
                edge_points.append(pts[-1])  # Store edge point for boundary layer
                pts.append(add_point(x3, y2))
                pts.append(add_point(x2, y1))
                pts.append(add_point(x1, y1))
                corner_points.append(pts[-1])  # Store corner point for refinement

                if i == N-1:
                    y4 = y1 + t
                    y3 = y4
                    Di, Db, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)

                    x0 = x1 - a
                    y0 = y1 + b

                    # ellipse points
                    pts.append(add_point(x0, y0))
                    pts.append(add_point(x0 + r, y0 + r))
                    pts.append(add_point(x0 + r, y0 + r -t))
                    pts.append(add_point(x0 + t, y0))


            pts.append(add_point(x1, y4))
            corner_points.append(pts[-1])  # Store corner point for refinement
            pts.append(add_point(x2, y4))
            pts.append(add_point(x3, y3))
            pts.append(add_point(xmax, y3))
            edge_points.append(pts[-1])  # Store edge point for boundary layer
    else:
        raise ValueError("Mw must be an even or odd integer.")
    
    # Add top points
    pts.append(add_point(xmax, ymax)) # top right
    pts.append(add_point(xmin, ymax)) # top left

    # Centers for arcs
    ell_center    = add_point(x1, y0)
    circ_center   = add_point(x0 + r, y0)
    circ_center2  = add_point(x0 + r, y0 + r - t/2)
    # -------------------------------

    # Build curves sequentially
    curve_tags = []
    for i in range(len(pts)-1 - 10):
        curve_tags.append(add_line(pts[i], pts[i+1]))
    
    # ellipse arcs
    curve_tags.append(add_ellipse_arc(pts[-11], ell_center, pts[-11], pts[-10]))
    curve_tags.append(add_circle_arc (pts[-10], circ_center,  pts[-9]))
    curve_tags.append(add_circle_arc (pts[-9], circ_center2, pts[-8]))
    curve_tags.append(add_circle_arc (pts[-8], circ_center,  pts[-7]))
    curve_tags.append(add_ellipse_arc(pts[-7], ell_center, pts[-7], pts[-6]))
    
    # last inlet lines
    for i in range(len(pts)-1 - 5, len(pts)-1):
        curve_tags.append(add_line(pts[i], pts[i+1]))
    # Add the last line to close the loop
    curve_tags.append(add_line(pts[-1], pts[0]))

    # Create surface
    loop_tag = gmsh.model.occ.addCurveLoop(curve_tags)
    surface  = gmsh.model.occ.addPlaneSurface([loop_tag])
    gmsh.model.occ.synchronize()

    return surface, curve_tags, edge_points, corner_points

# -----------------------------
# Boundary Layer
# -----------------------------

def boundaryLayerParameters(U_inf:float, nu:float=15.06e-6, x:float=0.3, y_plus:float=0.95, n_layers:int=20) -> tuple[float, float, float]:
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
    - y1: Height of the first cell [m]
    - expansion_ratio: Growth ratio between layers
    """
        
    if n_layers < 5:
        logging.warning("n_layers < 5, this may lead to inaccurate results.")

    # Calculate the first cell height (y1)
    Re_x = U_inf * x / nu                   # Reynolds number based on x 
    Cf = (2*math.log10(Re_x)-0.65)**(-2.3)  # Skin friction coefficient for turbulent boundary layer (empirical)
    u_tau = U_inf * math.sqrt(Cf / 2)       # Friction velocity
    yp = y_plus * nu / u_tau                # First cell height for desired y+
    y1 = 2*yp                               # Total height of the first cell 
    y1 = y1 / 7                             # Ajust from simulation results

    # Growth Ratio
    delta99 = 4.91*x/Re_x**0.5 if Re_x<5e5 else 0.38*x/Re_x**0.2 # Boundary layer thickness
    delta99 = delta99 * 1.00    # Adjust from simulation results
    r0 = 1.5                    # Initial guess for growth ratio
    
    def func(r:float) -> float:
        return  r**n_layers -r*(delta99/y1) + (delta99/y1-1)
    def fprime(r:float) -> float:
        return n_layers*r**(n_layers-1) - (delta99/y1)
    
    expansion_ratio = newton(func=func, x0=r0,fprime=fprime, tol=1e-4, maxiter=200)  # Solve for growth ratio using Newton's method
    # truncate expansion ration to two decimal places
    expansion_ratio = math.trunc(expansion_ratio*100)/100

    # Calculate turbulence frequency omega
    beta1 = 0.075
    omega = 6 * nu / (beta1 * (y1/2)**2)
    logging.info(f"turbulent frequency omega at wall: {omega:.2e} 1/s")

    return  y1, delta99, expansion_ratio

def apply_boundary_layer(curve_tags:list, edge_points:tuple[float, float], x:float, n_layers:int=20, y_plus:float=0.95,  U_inf:float=16) -> None:
    """
    Apply a structured boundary layer on selected curves and end points.

    Parameters:
    - curve_tags: List of curve tags for the boundary layer
    - edge_points: List of point tags for the fan points
    - x: Distance from leading edge [m] (for estimating Cf and delta99)
    - n_layers: Number of boundary layer layers
    - y_plus: Desired dimensionless first cell height
    - U_inf: Freestream velocity [m/s]
    """

    logging.info(f"edge_points : {edge_points}")
    logging.info(f"curve_tags  : {curve_tags}")
    
    y1, bl_thickness, expansion_ratio = boundaryLayerParameters(U_inf=U_inf, x=x, y_plus=y_plus, n_layers=n_layers)
    logging.info(f"Boundary layer parameters: y1={y1:.2e}, bl_thickness={bl_thickness:.2e}, expansion_ratio={expansion_ratio}, num_layers={n_layers}")

    bl = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(bl, "CurvesList", curve_tags)  # Curves for BL
    gmsh.model.mesh.field.setNumbers(bl, "PointsList", edge_points) # End points for BL
    gmsh.model.mesh.field.setNumber(bl, "Thickness", bl_thickness)  # Total thickness of boundary layer
    gmsh.model.mesh.field.setNumber(bl, "Size", y1)                 # Size of the first cell
    gmsh.model.mesh.field.setNumber(bl, "Ratio", expansion_ratio)   # Growth ratio between layers
    gmsh.model.mesh.field.setNumber(bl, "NbLayers", n_layers)       # Number of layers
    gmsh.model.mesh.field.setNumber(bl, "Quads", 1)                 # Use quadrilateral elements
    gmsh.model.mesh.field.setAsBoundaryLayer(bl)
    gmsh.model.occ.synchronize()

def apply_multiple_boundary_layers(all_edge_points: list, x: float, n_layers: int = 20, y_plus: float = 0.95, U_inf: float = 16) -> list:
    
    wall_curve_tags = []
    for i in range(0, len(all_edge_points), 2):
        edges = all_edge_points[i:i+2]
        
        if len(edges) < 2:
            logging.info(f"Skipping incomplete edge pair at index {i}: {edges}")
            continue
        
        curve_tags = list(range(edges[0], edges[1]))
        wall_curve_tags.extend(c for c in curve_tags)
        
        logging.info(f"Boundary layer {i//2+1}")
        apply_boundary_layer(
            curve_tags=curve_tags,
            edge_points=edges,
            x=x,
            n_layers=n_layers,
            y_plus=y_plus,
            U_inf=U_inf
        )
    gmsh.model.occ.synchronize() 
    return wall_curve_tags

# -----------------------------
# Sharp Corners
# ----------------------------
def refine_sharp_corner(point_tag:int, radius:float=1e-3, size:float=0.5e-3) -> None:
    """
    Refine sharp corners in the geometry by adding ball size field.
    """

    [X, Y, Z] = gmsh.model.getValue(0, point_tag, [])
    logging.info(f"Refining sharp corner at point {point_tag:2.0f} with coordinates ({X:.3f}, {Y:.3f}, {Z:.3f})")

    # Add ball size field for sharp corners
    ball = gmsh.model.mesh.field.add("Ball")
    gmsh.model.mesh.field.setNumber(ball, "Radius", radius)
    gmsh.model.mesh.field.setNumber(ball, "Thickness", 1e-2)
    gmsh.model.mesh.field.setNumber(ball, "VIn", size)
    gmsh.model.mesh.field.setNumber(ball, "VOut", 1e-2)
    gmsh.model.mesh.field.setNumber(ball, "XCenter", X)
    gmsh.model.mesh.field.setNumber(ball, "YCenter", Y)
    gmsh.model.mesh.field.setNumber(ball, "ZCenter", Z)

    # Set the ball size field as background mesh
    # gmsh.model.mesh.field.setAsBackgroundMesh(ball)
    gmsh.model.occ.synchronize()
    return ball

# -----------------------------
# Mesh Refinement
# -----------------------------
def refine_mesh(curve_tags:list, corner_points:list, xmax:float, ymin:float, Di:float, a:float, b:float, sizeThreshold:float=2e-3, sizeBox:float=2e-2) -> None:
    """
    Add distance + threshold + box fields and combine for smooth grading.
    """
    # Mesh zie from curvature
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 60)

    # Distance from walls
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "EdgesList", curve_tags)
    gmsh.model.mesh.field.setNumber(dist, "Sampling", 1e6)

    # Threshold based on distance
    thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thr, "InField", dist)
    gmsh.model.mesh.field.setNumber(thr, "SizeMin", sizeThreshold)
    gmsh.model.mesh.field.setNumber(thr, "SizeMax", 2)
    gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.03)
    gmsh.model.mesh.field.setNumber(thr, "DistMax", 5)

    # Box refinement near inlet corner
    XMin = xmax - L - a - 0.15
    YMax = ymin + Di/2 + b + 0.1
    logging.info(f"Box refinement near inlet corner: {XMin:.2f} < x < {xmax:.2f}, {ymin:.2f} < y < {YMax:.2f}")
    box = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(box, "VIn", sizeBox)
    gmsh.model.mesh.field.setNumber(box, "VOut", 1)
    gmsh.model.mesh.field.setNumber(box, "XMin", XMin)
    gmsh.model.mesh.field.setNumber(box, "XMax", xmax)
    gmsh.model.mesh.field.setNumber(box, "YMin", ymin)
    gmsh.model.mesh.field.setNumber(box, "YMax", YMax)
    gmsh.model.mesh.field.setNumber(box, "ZMin", -2)
    gmsh.model.mesh.field.setNumber(box, "ZMax",  2)
    gmsh.model.mesh.field.setNumber(box, "Thickness", 5)

    # Box refinement near inlet corner
    XMin = xmax - 0.01
    YMax = ymin + Di/2
    logging.info(f"Box refinement near inlets: {XMin:.2f} < x < {xmax:.2f}, {ymin:.2f} < y < {YMax:.2f}")
    box2 = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(box2, "VIn", sizeBox/2)
    gmsh.model.mesh.field.setNumber(box2, "VOut", 1)
    gmsh.model.mesh.field.setNumber(box2, "XMin", XMin)
    gmsh.model.mesh.field.setNumber(box2, "XMax", xmax)
    gmsh.model.mesh.field.setNumber(box2, "YMin", ymin)
    gmsh.model.mesh.field.setNumber(box2, "YMax", YMax)
    gmsh.model.mesh.field.setNumber(box2, "ZMin",  -0.01)
    gmsh.model.mesh.field.setNumber(box2, "ZMax",  1)
    gmsh.model.mesh.field.setNumber(box2, "Thickness", 0.00)

    fields = [thr, box, box2]

    # Refine sharp corners
    # for point_tag in corner_points[:-1]:
    #     # Refine sharp corners
    #     fields.append(refine_sharp_corner(point_tag=point_tag, radius=7e-3, size=0.2e-3))
    
    # Combine fields
    mfield = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(mfield, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(mfield)
    gmsh.model.occ.synchronize()

# -----------------------------
# Extrusion and Physical Groups
# -----------------------------
def extrude_and_group(surface:int, wall_curve_tags:list):
    """
    Extrude the 2D surface into a thin 3D volume by one element,
    then define physical groups for CFD boundaries and volume.
    """
    # Extrude surface into volume
    extruded = gmsh.model.occ.extrude([(2, surface)], 0, 0, 1.0,
                                       numElements=[1], recombine=True)
    gmsh.model.occ.synchronize()

    # get surface tags
    all_surface_dimTags = gmsh.model.getEntities(2)
    all_surface_tags = [tag for dim, tag in all_surface_dimTags]
    # print("all_surface_tags:", all_surface_tags)
    symmetryTag     = [all_surface_tags[1]]
    wallTag         = [c+1 for c in wall_curve_tags]
    inletTag        = [all_surface_tags[-4], all_surface_tags[-3], all_surface_tags[-2]]
    frontBackTag    = [all_surface_tags[-0], all_surface_tags[-1]]
    
    # Flatten all the tags we want to exclude
    excluded = set(symmetryTag + wallTag + inletTag + frontBackTag)
    # Compute outlet tags
    outletTag = [tag for tag in all_surface_tags if tag not in excluded]


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
    """
    Save the generated mesh to a file.
    """
    save_path.mkdir(parents=True, exist_ok=True)
    msh = save_path / f"{fname}.msh"
    gmsh.write(str(msh))

# -----------------------------
# Print Information
# -----------------------------    
def print_info(
        Mw:int, Mb:int, Kx:float, Ky:float, 
        r:float, 
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
    logging.info(f"          - Straight section len (mm)  : {L * 1e3:.1f}")
    logging.info(f"          - Bellmouth thickness (mm)   : {t * 1e3:.1f}")

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
         r:float=10e-3,
         xmin:float=0, ymin:float=0, xmax:float=25, ymax:float=25,
        nt:int=10, sd:str=None):
    """
    Main function to generate the mesh using GMSH.
    """

    # Define name and path for mesh file
    fname = f"ELL-{Mw}-{Mb}-{int(Kx*100)}-{int(Ky*100)}-{int(r*1e3)}-{int(L*1e3)}-{int(t*1e3)}"
    save_path = Path(sd) if sd else Path("outputs/meshes")
    # logging.info(f"Save path: {save_path}")
    

    # Initialize and add model
    gmsh.initialize()
    gmsh.model.add(fname)
    # num_threads = min(multiprocessing.cpu_count() // 2, nt)
    gmsh.option.setNumber("General.NumThreads", nt)     # Set number of threads for GMSH
    
    # Print info
    print_info(Mw=Mw, Mb=Mb, Kx=Kx, Ky=Ky, r=r,
                xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax,
                nt=nt, fname=fname)
    # Create geometry
    surf, _, edge_points, corner_points = create_geometry(Mw, Mb, Kx, Ky, r, xmin, ymin, xmax, ymax)

    # Apply boundary layer
    logging.info("Applying boundary layers...")
    wall_curve_tags = apply_multiple_boundary_layers(all_edge_points=edge_points, x=L, n_layers=25, y_plus=0.95, U_inf=16)
    
    # Refine mesh
    Di, _, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)
    logging.info(f"Inlet width Di:        {Di:.3f} m")
    logging.info(f"Inlet half width Di/2: {Di/2:.3f} m")
    logging.info("Refining mesh...")
    refine_mesh(curve_tags=wall_curve_tags, corner_points=corner_points, 
                xmax=xmax, ymin=ymin, Di=Di, a=a, b=b,
                sizeThreshold=2e-3, sizeBox=2e-2)

    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()

    # Generate 2D mesh
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.model.mesh.generate(2)

    # Extrude and define physical groups
    logging.info("Extruding and defining physical groups...")
    extrude_and_group(surface=surf,wall_curve_tags=wall_curve_tags)

    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()

    # Generate 3D mesh and save
    gmsh.model.mesh.generate(3)

    # Save mesh
    logging.info("Saving mesh...")
    save_mesh(fname, save_path=save_path)

    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()

    logging.info("Mesh generation complete. Finalizing GMSH...")
    gmsh.finalize()

    logging.info("Mesh finalized. Exiting...")

# -----------------------------
# Entry point for script execution
def run():
    args = parse_args()
    try:
        validate_args(args)
    except ValueError as e:
        logging.error(str(e))
        sys.exit(1)
    
    main(args.Mw, args.Mb, args.Kx, args.Ky, args.r,
         args.xmin, args.ymin, args.xmax, args.ymax,
         args.nt, args.sd)

if __name__ == "__main__":
    run()