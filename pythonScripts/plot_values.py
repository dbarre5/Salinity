# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

probeList = [10,20,90,150,220,295,360,380]
for data in range(len(probeList)):
    probeList[data] = (730-probeList[data])/100



#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

fileLocation = 'C:/Users/RDCHLDDB/Downloads/laminar_case_4_17/case.foam'

# create a new 'OpenFOAMReader'
casefoam = OpenFOAMReader(registrationName='case.foam', FileName=fileLocation)
casefoam.MeshRegions = ['internalMesh']
casefoam.CellArrays = ['U', 'alpha.sludge']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()


# Properties modified on animationScene1
animationScene1.AnimationTime = 1800.0


# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=casefoam)
contour1.ContourBy = ['POINTS', 'alpha.sludge']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

renderView1 = GetActiveViewOrCreate('RenderView')
# update the view to ensure updated data information
renderView1.Update()



# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=contour1)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [4.746389746665955, -1.5, 0.02536085993051529]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [4.746389746665955, -1.5, 0.02536085993051529]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 1.0, 0.0]



# update the view to ensure updated data information
renderView1.Update()

# save data
SaveData('C:/Users/RDCHLDDB/Downloads/laminar_case_4_17/contour.csv', proxy=slice1, PointDataArrays=['Normals', 'U', 'alpha.sludge', 'ddt0(alpha.sludge)', 'ddt0(rho,U)', 'ddtCorrDdt0(U)', 'multU', 'p', 'p_rgh'],
    CellDataArrays=['U', 'alpha.sludge', 'ddt0(alpha.sludge)', 'ddt0(rho,U)', 'ddtCorrDdt0(U)', 'multU', 'p', 'p_rgh'],
    FieldDataArrays=['CasePath'])


# get opacity transfer function/opacity map for 'alphasludge'
for dataPoint in range(len(probeList)):
    name = '_'+str(probeList[dataPoint])
    plotOverLineName = f'PlotOverLine{name}'

    # Create a *unique* PlotOverLine each loop
    pol = PlotOverLine(registrationName=plotOverLineName, Input=casefoam)

    # Set the line points
    x = probeList[dataPoint]
    pol.Point1 = [x, -1.5, 0.0]
    pol.Point2 = [x, -1.5, 0.1550000011920929]

    # Update
    renderView1.Update()


    # save data
    SaveData(name+".csv", proxy=pol, PointDataArrays=['U', 'alpha.sludge', 'arc_length', 'ddt0(alpha.sludge)', 'ddt0(rho,U)', 'ddtCorrDdt0(U)', 'multU', 'p', 'p_rgh', 'vtkValidPointMask'])


df = pd.read_csv("contour.csv")
x = df["Points:0"]*100
z = df["Points:2"]*100
plt.plot(x,z, marker="o", linestyle = "None", markersize=2)
velList = []
vertLocList = []
horLocList = []

for dataPoint in range(len(probeList)):
    name = '_'+str(probeList[dataPoint])
    filename = name+".csv"
    #Read data for vertical lines
    df = pd.read_csv(filename)
    print(df)
    velList.append(df["U:0"])
    vertLocList.append(df["Points:2"]*100)
    print(vertLocList[0])
    horLocList.append(df["U:0"]*500+df["Points:0"]*100)
    print(horLocList[0])



    figure(figsize=(20, 4), dpi=300)

    plt.plot(horLocList[dataPoint],vertLocList[dataPoint], marker="o", linestyle="None", markersize=2)
plt.ylim(0,16)
plt.xlim(0,760)
plt.savefig("testPlot.png")
