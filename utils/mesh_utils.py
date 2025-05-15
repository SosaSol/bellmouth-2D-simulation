from utils.config import *
import gmsh
from scipy.optimize import newton
from pathlib import Path
import math
import logging

# -----------------------------
# Boundary Layer
# -----------------------------

def boundaryLayerParameters(U_inf:float, nu:float=15.06e-6, x:float=0.3, y_plus:float=0.95, n_layers:int=20, yp_factor:float=1) -> tuple[float, float, float]:
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
    - delta99: Boundary layer thickness [m]
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
    y1 = y1 * yp_factor                     # Ajust from simulation results

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

def apply_boundary_layer(curve_tags:list, edge_points:tuple[float, float], x:float, n_layers:int=20, y_plus:float=0.95,  U_inf:float=16, yp_factor:float=1) -> None:
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
    
    y1, bl_thickness, expansion_ratio = boundaryLayerParameters(U_inf=U_inf, x=x, y_plus=y_plus, n_layers=n_layers, yp_factor=yp_factor)
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

def apply_multiple_boundary_layers(all_edge_points: list, x: float, n_layers: int = 20, y_plus: float = 0.95, U_inf: float = 16, yp_factor:float=1) -> list:
    """
    Apply multiple boundary layers to the geometry based on edge points.
    
    Parameters:
    - all_edge_points (list): List of edge points for the boundary layers
    - x (float): Distance from leading edge [m] (for estimating Cf and delta99)
    - n_layers (int): Number of boundary layer layers
    - y_plus (float): Desired dimensionless first cell height
    - U_inf (float): Freestream velocity [m/s]
   
    Returns:
    - wall_curve_tags: List of curve tags for the wall boundary layers
    """

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
            U_inf=U_inf,
            yp_factor=yp_factor
        )
    gmsh.model.occ.synchronize() 
    return wall_curve_tags

# -----------------------------
# Sharp Corners
# ----------------------------
def refine_sharp_corner(point_tag:int, radius:float=1e-3, size:float=0.5e-3) -> None:
    """
    Refine sharp corners in the geometry by adding ball size field.
    NEVER USED
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
def refine_mesh(curve_tags:list, xmax:float, ymin:float, Di:float, a:float, b:float, sizeThreshold:float=2e-3, sizeBox:float=2e-2) -> None:
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
    gmsh.model.mesh.field.setNumber(thr, "SizeMax", 1)
    gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.03)
    gmsh.model.mesh.field.setNumber(thr, "DistMax", 0.03)

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
    gmsh.model.mesh.field.setNumber(box, "ZMin", -0.01)
    gmsh.model.mesh.field.setNumber(box, "ZMax", 1)
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
    gmsh.model.occ.extrude([(2, surface)], 0, 0, 1.0, numElements=[1], recombine=True)
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
def save_mesh(fname:str , save_path:Path):
    """
    Save the generated mesh to a file.
    """
    save_path.mkdir(parents=True, exist_ok=True)
    msh = save_path / f"{fname}.msh"
    gmsh.write(str(msh))