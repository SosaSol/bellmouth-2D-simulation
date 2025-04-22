import gmsh
import sys
import multiprocessing
import argparse
import os

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
    return parser.parse_args()

# -----------------------------
# Geometry Helper Functions
# -----------------------------
def compute_geometry_parameters(Mw, Mb, Kx, Ky):
    """
    Compute key geometry parameters based on input multipliers.
    Returns:
        Di: inlet width
        Db: bellmouth reference width
        a: semi-major axis (half bellmouth width)
        b: semi-minor axis (half bellmouth height)
    """
    Di = 242.09e-3 * Mw + 3.5455e-3 # inlet width
    Db = 242.09e-3 * Mb + 3.5455e-3 # bellmouth reference width
    # Ensure a >= b for ellipse definition
    a_raw = Db / 2
    b_raw = Ky * Db / 2
    a = max(a_raw, b_raw)
    b = min(a_raw, b_raw)
    return Di, Db, a, b


def add_point(x, y, z=0):
    """Add a point in the OCC geometry kernel."""
    return gmsh.model.occ.addPoint(x, y, z)


def add_line(p1, p2):
    """Add a straight line between two points."""
    return gmsh.model.occ.addLine(p1, p2)


def add_circle_arc(p1, center, p2):
    """Add a circular arc from p1 to p2 around center."""
    return gmsh.model.occ.addCircleArc(p1, center, p2)


def add_ellipse_arc(p1, center, major, p2):
    """Add an elliptical arc from p1 to p2 using major axis vector."""
    return gmsh.model.occ.addEllipseArc(p1, center, major, p2)

# -----------------------------
# Build 2D Geometry
# -----------------------------
def create_geometry(Mw, Mb, Kx, Ky, r, t, L, xmin, ymin, xmax, ymax):
    """
    Build the parametric 2D bellmouth inlet geometry.
    Returns:
        surface: the tag of the created surface
        curve_tags: list of all curve tags in order
        edge_points: two point tags for boundary layer 'fan'
    """
    print("Info    : Building geometry...")

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
def apply_boundary_layer(curve_tags, edge_points, bl_thickness=3e-3):
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
def refine_mesh(curve_tags, xmin, xmax, ymin, ymax):
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
def extrude_and_group(surface):
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
def save_mesh(fname, save_path):
    """
    Save the final mesh with a filename encoding key parameters
    in the specified directory.
    """
    # Ensure the output directory exists
    os.makedirs(save_path, exist_ok=True)

    # Construct full path
    msh_name = os.path.join(save_path, f"{fname}.msh")
    
    # Save the mesh
    gmsh.write(msh_name)
    print(f"Info    : Saving mesh to {save_path}")
    print(f"Mesh saved to {msh_name}")


# -----------------------------
# Print Information
# -----------------------------    
def print_info(Mw, Mb, Kx, Ky, r, t, L, xmin, ymin, xmax, ymax, nt, fname):
    """
    Print the parameters used for mesh generation.
    """
    print("Info    : Geometry parameters:")
    print(f"          - Mesh tag suffix            : {fname}")
    print(f"          - WindShaper modules         : {Mw}")
    print(f"          - Bellmouth reference modules: {Mb}")
    print(f"          - Ellipse scaling (Kx, Ky)   : {Kx:.2f}, {Ky:.2f}")
    print(f"          - Radius of curvature (mm)   : {r * 1e3:.1f}")
    print(f"          - Wall thickness (mm)        : {t * 1e3:.1f}")
    print(f"          - Straight section len (mm)  : {L * 1e3:.1f}")

    print(f"Info    : Domain bounds -> X: [{xmin}, {xmax}], Y: [{ymin}, {ymax}]")
    print(f"Info    : Using {nt} threads for GMSH.")

# -----------------------------
# Main Execution
# -----------------------------
def main(Mw=12, Mb=12, Kx=0.33, Ky=0.33, r=10e-3, t=5e-3, L=0.3, xmin=0, ymin=0, xmax=25, ymax=25, nt=10):
    """
    Main function to generate the mesh using GMSH.
    """

    # Define name and path for mesh file
    fname = f"ELL-{Mw}-{Mb}-{int(Kx*100)}-{int(Ky*100)}-{int(r*1e3)}-{int(t*1e3)}"
    save_path = r"\\wsl.localhost\Ubuntu\home\solim\OpenFOAM\solim-v2412\run\bellmouth_2D\{}\constant\triSurface".format(fname)
    # print(f"Info    : Save path: {save_path}")

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
if __name__ == "__main__":
    args = parse_args()
    main(args.Mw, args.Mb, args.Kx, args.Ky, args.r, args.t,
         args.L, args.xmin, args.ymin, args.xmax, args.ymax,
         args.nt)