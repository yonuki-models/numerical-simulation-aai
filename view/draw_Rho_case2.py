# state file generated using paraview version 5.9.1
import os

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

with open('./log/log_Rho.txt', 'a') as f:
    f.write('test A\n')

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
# renderView1.ViewSize = [2296, 880]
renderView1.ViewSize = [1800, 900]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.0, 0.0, -100.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2000.0, -3800.0, 2000.0]
renderView1.CameraFocalPoint = [100.0, 0.0, -400.0]
renderView1.CameraViewUp = [
    -0.17443747914703292, 0.34887495829406584, 0.9207919576886748]
renderView1.CameraViewAngle = 50.0
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 2400
renderView1.CameraParallelProjection = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
# layout1.SetSize(2296, 880)
layout1.SetSize(1800, 900)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Partitioned Structured Grid Reader'
experiment_name = "exRo9N10e01_LES"
uVWT000 = XMLPartitionedStructuredGridReader(
    registrationName='UVWT000*',
    FileName=['/work/03/gs53/c24070/spectral_model_distortion_SSL/data/'
              + experiment_name
              + '/vtk/UVWT{:0>6}.pvts'.format(i) for i in range(0, 280)]
)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

uVWT000.PointArrayStatus = ['T', 'U']
uVWT000.TimeArray = 'None'

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=uVWT000)
calculator1.CoordinateResults = 1
calculator1.Function = '-coordsX*iHat-coordsY*jHat+5*coordsZ*kHat'

# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=calculator1)
calculator3.CoordinateResults = 1
calculator3.Function = 'coordsX*iHat+coordsY*jHat+(coordsZ-1000)*kHat'

# create a new 'Calculator'
calculator4 = Calculator(registrationName='Calculator4', Input=calculator3)
calculator4.ResultArrayName = 'Rho'
calculator4.Function = 'T+0.0002*coordsZ'

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)
calculator2.ResultArrayName = 'Rho'
calculator2.Function = 'T+0.0002*coordsZ'

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=calculator2)
contour1.ContourBy = ['POINTS', 'Rho']
contour1.ComputeNormals = 0
contour1.GenerateTriangles = 0
contour1.Isosurfaces = [1e-08]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=calculator2)
threshold1.Scalars = ['POINTS', 'Rho']
threshold1.ThresholdRange = [-0.2, 0.0]

# create a new 'Threshold'
threshold2 = Threshold(registrationName='Threshold2', Input=calculator4)
threshold2.Scalars = ['POINTS', 'Rho']
threshold2.ThresholdRange = [-0.2, -0.0]

# create a new 'Ellipse'
# ellipse1 = Ellipse(registrationName='Ellipse1')
# ellipse1.Center = [0.0, 0.0, -1000.0]
# ellipse1.MajorRadiusVector = [2234.3573854970023, 0.0, 0.0]
# ellipse1.Ratio = 0.4

# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=calculator4)
contour2.ContourBy = ['POINTS', 'Rho']
contour2.ComputeNormals = 0
contour2.GenerateTriangles = 0
contour2.Isosurfaces = [-0.2000001]
contour2.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Rho'
rhoLUT = GetColorTransferFunction('Rho')
rhoPWF = GetOpacityTransferFunction('Rho')
rhoLUT.ApplyPreset('Rainbow Desaturated', True)
contour1Display.RescaleTransferFunctionToDataRange(False, True)
rhoLUT.RescaleTransferFunction(-0.2, 0.0)
rhoPWF.RescaleTransferFunction(-0.2, 0.0)

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'Rho']
contour1Display.LookupTable = rhoLUT
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'None'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'T'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'U'
contour1Display.ScaleFactor = 315.9858642578125
contour1Display.SelectScaleArray = 'None'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'None'
contour1Display.GaussianRadius = 15.799293212890625
contour1Display.SetScaleArray = ['POINTS', 'T']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'T']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# contour1Display.ScaleTransferFunction.Points = [
#     -0.016501646488904953, 0.0, 0.5, 0.0,
#     0.02443234995007515, 1.0, 0.5, 0.0
# ]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# contour1Display.OpacityTransferFunction.Points = [
#     -0.016501646488904953, 0.0, 0.5, 0.0,
#     0.02443234995007515, 1.0, 0.5, 0.0
# ]

# show data from threshold1
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# get opacity transfer function/opacity map for 'Rho'
# rhoPWF = GetOpacityTransferFunction('Rho')
# rhoPWF.Points = [-0.061585888266563416, 0.0, 0.5, 0.0, 1.0001818928628836e-08, 1.0, 0.5, 0.0]
# rhoPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['POINTS', 'Rho']
threshold1Display.LookupTable = rhoLUT
threshold1Display.SelectTCoordArray = 'None'
threshold1Display.SelectNormalArray = 'None'
threshold1Display.SelectTangentArray = 'None'
threshold1Display.OSPRayScaleArray = 'Rho'
threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold1Display.SelectOrientationVectors = 'U'
threshold1Display.ScaleFactor = 315.9858642578125
threshold1Display.SelectScaleArray = 'Rho'
threshold1Display.GlyphType = 'Arrow'
threshold1Display.GlyphTableIndexArray = 'Rho'
threshold1Display.GaussianRadius = 15.799293212890625
threshold1Display.SetScaleArray = ['POINTS', 'Rho']
threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold1Display.OpacityArray = ['POINTS', 'Rho']
threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
threshold1Display.PolarAxes = 'PolarAxesRepresentation'
threshold1Display.ScalarOpacityFunction = rhoPWF
threshold1Display.ScalarOpacityUnitDistance = 8.973830108657838
threshold1Display.OpacityArrayName = ['POINTS', 'Rho']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# threshold1Display.ScaleTransferFunction.Points = [-0.05412979291677475, 0.0, 0.5, 0.0, -6.064772573033395e-10, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# threshold1Display.OpacityTransferFunction.Points = [-0.05412979291677475, 0.0, 0.5, 0.0, -6.064772573033395e-10, 1.0, 0.5, 0.0]

# show data from threshold2
threshold2Display = Show(threshold2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'
threshold2Display.ColorArrayName = ['POINTS', 'Rho']
threshold2Display.LookupTable = rhoLUT
threshold2Display.SelectTCoordArray = 'None'
threshold2Display.SelectNormalArray = 'None'
threshold2Display.SelectTangentArray = 'None'
threshold2Display.OSPRayScaleArray = 'Rho'
threshold2Display.OSPRayScaleFunction = 'PiecewiseFunction'
threshold2Display.SelectOrientationVectors = 'U'
threshold2Display.ScaleFactor = 315.9858642578125
threshold2Display.SelectScaleArray = 'Rho'
threshold2Display.GlyphType = 'Arrow'
threshold2Display.GlyphTableIndexArray = 'Rho'
threshold2Display.GaussianRadius = 15.799293212890625
threshold2Display.SetScaleArray = ['POINTS', 'Rho']
threshold2Display.ScaleTransferFunction = 'PiecewiseFunction'
threshold2Display.OpacityArray = ['POINTS', 'Rho']
threshold2Display.OpacityTransferFunction = 'PiecewiseFunction'
threshold2Display.DataAxesGrid = 'GridAxesRepresentation'
threshold2Display.PolarAxes = 'PolarAxesRepresentation'
threshold2Display.ScalarOpacityFunction = rhoPWF
threshold2Display.ScalarOpacityUnitDistance = 8.97526449677581
threshold2Display.OpacityArrayName = ['POINTS', 'Rho']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# threshold2Display.ScaleTransferFunction.Points = [-0.0599999995562629, 0.0, 0.5, 0.0, -0.009041980579495436, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# threshold2Display.OpacityTransferFunction.Points = [-0.0599999995562629, 0.0, 0.5, 0.0, -0.009041980579495436, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'Rho']
contour2Display.LookupTable = rhoLUT
contour2Display.SelectTCoordArray = 'None'
contour2Display.SelectNormalArray = 'None'
contour2Display.SelectTangentArray = 'None'
contour2Display.OSPRayScaleArray = 'Rho'
contour2Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour2Display.SelectOrientationVectors = 'U'
contour2Display.ScaleFactor = 315.9858642578125
contour2Display.SelectScaleArray = 'Rho'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'Rho'
contour2Display.GaussianRadius = 15.799293212890625
contour2Display.SetScaleArray = ['POINTS', 'Rho']
contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
contour2Display.OpacityArray = ['POINTS', 'Rho']
contour2Display.OpacityTransferFunction = 'PiecewiseFunction'
contour2Display.DataAxesGrid = 'GridAxesRepresentation'
contour2Display.PolarAxes = 'PolarAxesRepresentation'

# create a new 'Ellipse'
ellipse1 = Ellipse(registrationName='Ellipse1')
ellipse1.Center = [0.0, 0.0, -1000.0]
ellipse1.MajorRadiusVector = [3422.4, 0.0, 0.0]
                              # 200 * 24.20 / 2**0.5
ellipse1.Ratio = 0.9

# show data from ellipse1
ellipse1Display = Show(ellipse1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
ellipse1Display.Representation = 'Surface'
ellipse1Display.AmbientColor = [0.0, 0.0, 0.0]
ellipse1Display.ColorArrayName = [None, '']
ellipse1Display.DiffuseColor = [0.0, 0.0, 0.0]
ellipse1Display.LineWidth = 2.0
ellipse1Display.SelectTCoordArray = 'Texture Coordinates'
ellipse1Display.SelectNormalArray = 'None'
ellipse1Display.SelectTangentArray = 'None'
ellipse1Display.OSPRayScaleArray = 'Texture Coordinates'
ellipse1Display.OSPRayScaleFunction = 'PiecewiseFunction'
ellipse1Display.SelectOrientationVectors = 'None'
ellipse1Display.ScaleFactor = 446.871484375
ellipse1Display.SelectScaleArray = 'None'
ellipse1Display.GlyphType = 'Arrow'
ellipse1Display.GlyphTableIndexArray = 'None'
ellipse1Display.GaussianRadius = 22.34357421875
ellipse1Display.SetScaleArray = ['POINTS', 'Texture Coordinates']
ellipse1Display.ScaleTransferFunction = 'PiecewiseFunction'
ellipse1Display.OpacityArray = ['POINTS', 'Texture Coordinates']
ellipse1Display.OpacityTransferFunction = 'PiecewiseFunction'
ellipse1Display.DataAxesGrid = 'GridAxesRepresentation'
ellipse1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
ellipse1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9900000095367432, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
ellipse1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9900000095367432, 1.0, 0.5, 0.0]


# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# contour2Display.ScaleTransferFunction.Points = [-0.061585888266563416, 0.0, 0.5, 0.0, -0.061578258872032166, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# contour2Display.OpacityTransferFunction.Points = [-0.061585888266563416, 0.0, 0.5, 0.0, -0.061578258872032166, 1.0, 0.5, 0.0]

# show data from ellipse1
# ellipse1Display = Show(ellipse1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
# ellipse1Display.Representation = 'Surface'
# ellipse1Display.AmbientColor = [0.0, 0.0, 0.0]
# ellipse1Display.ColorArrayName = [None, '']
# ellipse1Display.DiffuseColor = [0.0, 0.0, 0.0]
# ellipse1Display.LineWidth = 2.0
# ellipse1Display.SelectTCoordArray = 'Texture Coordinates'
# ellipse1Display.SelectNormalArray = 'None'
# ellipse1Display.SelectTangentArray = 'None'
# ellipse1Display.OSPRayScaleArray = 'Texture Coordinates'
# ellipse1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# ellipse1Display.SelectOrientationVectors = 'None'
# ellipse1Display.ScaleFactor = 446.871484375
# ellipse1Display.SelectScaleArray = 'None'
# ellipse1Display.GlyphType = 'Arrow'
# ellipse1Display.GlyphTableIndexArray = 'None'
# ellipse1Display.GaussianRadius = 22.34357421875
# ellipse1Display.SetScaleArray = ['POINTS', 'Texture Coordinates']
# ellipse1Display.ScaleTransferFunction = 'PiecewiseFunction'
# ellipse1Display.OpacityArray = ['POINTS', 'Texture Coordinates']
# ellipse1Display.OpacityTransferFunction = 'PiecewiseFunction'
# ellipse1Display.DataAxesGrid = 'GridAxesRepresentation'
# ellipse1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# ellipse1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9900000095367432, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# ellipse1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.9900000095367432, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for rhoLUT in view renderView1
rhoLUTColorBar = GetScalarBar(rhoLUT, renderView1)
rhoLUTColorBar.Title = 'Rho'
rhoLUTColorBar.ComponentTitle = ''

# set color bar visibility
# rhoLUTColorBar.Visibility = 1
rhoLUTColorBar.Visibility = 0

# show color legend
# contour1Display.SetScalarBarVisibility(renderView1, False)

# show color legend
# threshold1Display.SetScalarBarVisibility(renderView1, False)

# show color legend
# threshold2Display.SetScalarBarVisibility(renderView1, False)

# show color legend
# contour2Display.SetScalarBarVisibility(renderView1, False)

# rhoLUTColorBar.Visibility = 0

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(None)
# ----------------------------------------------------------------

for i in range(0, 280):
    animationScene1.AnimationTime = i * 1.e0
    anime_directory = "../anime/" + experiment_name
    if not os.path.isdir(anime_directory):
        os.makedirs(anime_directory)
    SaveScreenshot(
        anime_directory + '/Rho{:0>6}.png'.format(i),
        renderView1
    )
    with open('./log/log_Rho.txt', 'a') as f:
        f.write('Rho ' + str(i) + '\n')

if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
    exit()
