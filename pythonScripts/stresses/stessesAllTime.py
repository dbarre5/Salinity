from paraview.simple import *
import os

try:
    base_path = os.path.dirname(os.path.abspath(__file__))
except NameError:
    base_path = os.getcwd()

fileName = os.path.join(base_path, "case.foam")
casefoam = OpenFOAMReader(FileName=fileName)

# Make sure we can access all timesteps
casefoam.UpdatePipeline()
timesteps = casefoam.TimestepValues
print("Available timesteps:", timesteps)

# Example of processing a derived filter (your CellCenters filter)
gradient1 = Gradient(Input=casefoam)
gradient1.ScalarArray = ['CELLS', 'U']

extractComponent1 = ExtractComponent(Input=gradient1)
extractComponent1.InputArray = ['POINTS', 'U']
extractComponent1.OutputArrayName = 'Ux'

contour2 = Contour(Input=extractComponent1)
contour2.ContourBy = ['POINTS', 'Ux']
contour2.Isosurfaces = [0.0]


meshQuality1 = MeshQuality(Input=contour2)
meshQuality1.TriangleQualityMeasure = 'Area'
meshQuality1.QuadQualityMeasure = 'Area'

cellCenters1 = CellCenters(Input=meshQuality1)
cellCenters1.VertexCells = 1

# Loop through all timesteps and save a CSV for each
for t in timesteps:
    print(f"Processing timestep: {t}")
    
    # Set the timestep
    casefoam.UpdatePipeline(time=t)
    
    # Update filters
    gradient1.UpdatePipeline(time=t)
    extractComponent1.UpdatePipeline(time=t)
    contour2.UpdatePipeline(time=t)
    meshQuality1.UpdatePipeline(time=t)
    cellCenters1.UpdatePipeline(time=t)
    
    # Save CSV
    csv_path = os.path.join(base_path, f"u_t{t:.6f}.csv")
    SaveData(csv_path, proxy=cellCenters1)
    print(f"Saved CSV for timestep {t} at {csv_path}")



# Example of processing a derived filter (your CellCenters filter)
gradient2 = Gradient(Input=casefoam)
gradient2.ScalarArray = ['CELLS', 'p_rgh']

extractComponent2 = ExtractComponent(Input=gradient2)
extractComponent2.InputArray = ['POINTS', 'U']
extractComponent2.OutputArrayName = 'Ux'
###
contour3 = Contour(Input=extractComponent2)
contour3.ContourBy = ['POINTS', 'Ux']
contour3.Isosurfaces = [0.0]

generateSurfaceNormals1 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=contour3)
# Properties modified on generateSurfaceNormals1
generateSurfaceNormals1.ComputeCellNormals = 1

meshQuality2 = MeshQuality(Input=generateSurfaceNormals1)
meshQuality2.TriangleQualityMeasure = 'Area'
meshQuality2.QuadQualityMeasure = 'Area'

cellCenters2 = CellCenters(Input=meshQuality2)
cellCenters2.VertexCells = 1

# Loop through all timesteps and save a CSV for each
for t in timesteps:
    print(f"Processing timestep: {t}")
    
    # Set the timestep
    casefoam.UpdatePipeline(time=t)
    
    # Update filters
    gradient2.UpdatePipeline(time=t)
    extractComponent2.UpdatePipeline(time=t)
    contour3.UpdatePipeline(time=t)
    meshQuality2.UpdatePipeline(time=t)
    cellCenters2.UpdatePipeline(time=t)
    
    # Save CSV
    csv_path = os.path.join(base_path, f"p_t{t:.6f}.csv")
    SaveData(csv_path, proxy=cellCenters2)
    print(f"Saved CSV for timestep {t} at {csv_path}")
