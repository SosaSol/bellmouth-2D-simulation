#!/usr/bin/env python3
# ==============================================================================
# mesh_advanced_with_bellmouth.py
#
# 2D Geometry Generator for WindShaper "Aeropack" Inlet
# -----------------------------------------------------
# This script generates the 2D geometry for the WindShaper Aeropack,
# including only the bellmouth.
# The geometry is fully parametric and suitable for advanced mesh generation
# with boundary layers and refinement zones.
#
# Author: Solim Rovera
# Date:   14-05-2025
# ==============================================================================

# import libraries
import sys
import argparse
import gmsh # GMSH library

# Import config
from utils.config import *
from utils.geometry_utils import *
from utils.mesh_utils import *

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
    
    parser.add_argument("--t", type=float, default=10e-3, help="Bellmouth thickness (in meters)")
    parser.add_argument("--r", type=float, default=100e-3, help="Radius of curvature (in meters)")

    parser.add_argument("--xmin", type=float, default=0, help="X min (in meters)")
    parser.add_argument("--ymin", type=float, default=0, help="Y min (in meters)")
    parser.add_argument("--xmax", type=float, default=25, help="X max (in meters)")
    parser.add_argument("--ymax", type=float, default=25, help="Y max (in meters)")

    parser.add_argument("--nt", type=int, default=10, help="Number of threads for GMSH")
    parser.add_argument("--sd", type=str, default=None, help="Save directory")
    return parser.parse_args()

# -----------------------------
# Build 2D Geometry
# -----------------------------
def create_geometry(
        Mw:int, Mb:int, Kx: float, Ky:float, 
        t:float, r:float, 
        xmin:float, ymin:float, xmax:float, ymax:float
        ) -> tuple[int, list[int], list[int]]:
    """
    Build the parametric 2D bellmouth inlet geometry.
    Parameters:
        Mw (int): Number of WindShaper modules
        Mb (int): Bellmouth reference modules
        Kx (float): X-axis ellipse scale
        Ky (float): Y-axis ellipse scale
        t (float): Bellmouth thickness
        r (float): Radius of curvature
        xmin (float): Minimum x-coordinate
        ymin (float): Minimum y-coordinate
        xmax (float): Maximum x-coordinate
        ymax (float): Maximum y-coordinate
    Returns:
        surface (int): Surface tag of the created geometry
        curve_tags (list[int]): List of curve tags for the geometry
        edge_points (list[int]): List of edge point tags for boundary layer
        corner_points (list[int]): List of corner point tags for refinement
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

    # define the end thickness t_end
    if Mw <= 3:
        t_end = t_22
    elif Mw <= 12:
        t_end = t_44
    else:
        # raise an error saying that the thickness is not defined for this number of modules
        raise ValueError(f"t_end not defined for {Mw} modules.")

    # If Mw odd
    if Mw % 2 == 1:
        pts.append(add_point(xmax, ymin))            # p2         
        N = (Mw + 1) // 2
        for i in range(N):
            y1 = H1 * (i+1/2) + GAP * i
            y2 = y1 + DELTA_H
            y4 = y1 + GAP
            y3 = y4 - DELTA_H

            pts.append(add_point(xmax, y2))
            edge_points.append(pts[-1])  # Store edge point for boundary layer
            pts.append(add_point(x3, y2))
            pts.append(add_point(x2, y1))
            pts.append(add_point(x1, y1))
            corner_points.append(pts[-1])  # Store corner point for refinement

            if i == N-1:
                y4 = y1 + t_end
                y3 = y4
                Di, Db, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)

                x0 = x1 - a
                y0 = y1 + b

                # ellipse points
                pts.append(add_point(x0, y0))
                pts.append(add_point(x0 + r, y0 + r))
                pts.append(add_point(x0 + r, y0 + r -t))
                pts.append(add_point(x0 + t, y0))
                pts.append(add_point(x1-0.05, y1 + t)) # -0.05 because the Boundary Layer breaks


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
            y1 = H1 * i + GAP * (i-1/2)
            y2 = y1 + DELTA_H
            y4 = y1 + GAP
            y3 = y4 - DELTA_H

            if i == 0:
                pts.append(add_point(x1, ymin)) # p2
                edge_points.append(pts[-1])     # Store edge point for boundary layer
                y4 = GAP/2
                y3 = y4 - DELTA_H
            elif i!= 0:
                pts.append(add_point(xmax, y2))
                edge_points.append(pts[-1])     # Store edge point for boundary layer
                pts.append(add_point(x3, y2))
                pts.append(add_point(x2, y1))
                pts.append(add_point(x1, y1))
                corner_points.append(pts[-1])   # Store corner point for refinement

                if i == N-1:
                    y4 = y1 + t_end
                    y3 = y4
                    Di, Db, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)

                    x0 = x1 - a
                    y0 = y1 + b

                    # ellipse points
                    pts.append(add_point(x0, y0))
                    pts.append(add_point(x0 + r, y0 + r))
                    pts.append(add_point(x0 + r, y0 + r -t))
                    pts.append(add_point(x0 + t, y0))
                    pts.append(add_point(x1-0.05, y1 + t)) # -0.05 because the Boundary Layer breaks


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
    for i in range(len(pts)-1 - 11):
        curve_tags.append(add_line(pts[i], pts[i+1]))
    
    # bellmouth ellipse arcs
    curve_tags.append(add_ellipse_arc(pts[-12], ell_center, pts[-12], pts[-11]))
    curve_tags.append(add_circle_arc (pts[-11], circ_center,  pts[-10]))
    curve_tags.append(add_circle_arc (pts[-10], circ_center2, pts[-9]))
    curve_tags.append(add_circle_arc (pts[-9], circ_center,  pts[-8]))
    curve_tags.append(add_ellipse_arc(pts[-8], ell_center, pts[-8], pts[-7]))
    
    # last inlet lines
    for i in range(len(pts)-1 - 6, len(pts)-1):
        curve_tags.append(add_line(pts[i], pts[i+1]))
    # Add the last line to close the loop
    curve_tags.append(add_line(pts[-1], pts[0]))

    # Create surface
    loop_tag = gmsh.model.occ.addCurveLoop(curve_tags)
    surface  = gmsh.model.occ.addPlaneSurface([loop_tag])
    gmsh.model.occ.synchronize()

    return surface, curve_tags, edge_points, corner_points

# -----------------------------
# Print Information
# -----------------------------    
def print_info(
        Mw:int, Mb:int, Kx:float, Ky:float, 
        t:float, r:float, 
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
    logging.info(f"          - Bellmouth thickness (mm)   : {t * 1e3:.1f}")
    logging.info(f"          - Radius of curvature (mm)   : {r * 1e3:.1f}")

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
        raise ValueError("Mw (Number of WindShaper modules) must be greater than 0.")
    if args.Mb <= 0:
        raise ValueError("Mb (Bellmouth reference modules) must be greater than 0.")
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
def main(Mw:int, Mb:int, Kx:float, Ky:float, 
         t:float, r:float, 
         xmin:float, ymin:float, xmax:float, ymax:float,
         nt:int, sd:str):
    """
    Main function to generate the mesh using GMSH.
    """

    # Define name and path for mesh file
    fname = f"ELL-{Mw}-{Mb}-{int(Kx*100)}-{int(Ky*100)}-{int(t*1e3)}-{int(r*1e3)}"
    save_path = Path(sd) if sd else Path("outputs/meshes")
    # logging.info(f"Save path: {save_path}")
    

    # Initialize and add model
    gmsh.initialize()
    gmsh.model.add(fname)
    # num_threads = min(multiprocessing.cpu_count() // 2, nt)
    gmsh.option.setNumber("General.NumThreads", nt)     # Set number of threads for GMSH
    
    # Print info
    print_info(Mw=Mw, Mb=Mb, Kx=Kx, Ky=Ky, t=t, r=r,
                xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax,
                nt=nt, fname=fname)
    # Create geometry
    surf, _, edge_points, _ = create_geometry(Mw, Mb, Kx, Ky, t, r, xmin, ymin, xmax, ymax)

    # Apply boundary layer
    logging.info("Applying boundary layers...")
    wall_curve_tags = apply_multiple_boundary_layers(all_edge_points=edge_points, x=L, n_layers=21,
                                                      y_plus=0.95, U_inf=16, yp_factor=0.3)
    
    # Refine mesh
    Di, Db, a, b = compute_geometry_parameters(Mw, Mb, Kx, Ky)
    logging.info(f"Inlet width Di:        {Di:.3f} m")
    logging.info(f"Inlet half width Di/2: {Di/2:.3f} m")
    logging.info(f"Bellmouth semi-m axis: {Db:.3f} m")
    logging.info("Refining mesh...")
    refine_mesh(curve_tags=wall_curve_tags, 
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
    
    main(args.Mw, args.Mb, args.Kx, args.Ky, args.t, args.r,
         args.xmin, args.ymin, args.xmax, args.ymax,
         args.nt, args.sd)

if __name__ == "__main__":
    run()