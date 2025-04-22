# 


# 1. Create the mesh
cartesian2DMesh # 2D mesh
# flattenMesh # Flatten the mesh for better visualization

# foamMeshToFluent

# 2. Check the mesh
checkMesh
# foamToVTK -pointSet nonAlignedEdges
# for point sets, use the Glyph filter in ParaView, with spheres and adjust the size

# foamToVTK -faceSet nonOrthoFaces -time 0


# 2. Run the simulation
simpleFoam

# 3. Check the y+ condition
postProcess -func yPlus

# 4. Post-process the results
postProcess -func "patchIntegrate(name=outlet,U)" 
postProcess -func "patchAverage(name=inlet,p)"
postProcess -func "patchAverage(name=outlet,p)"


# 3. Post-process the results
touch foam.foam
