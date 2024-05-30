# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023 replay file
# Internal Version: 2022_09_28-20.11.55 183150
# Run by cnober on Thu May 30 14:44:45 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=375.819793701172, 
    height=212.760009765625)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='UT_EBT_polyN.odb')
#: Model: C:/temp/polyN_optimization/execution/UT_EBT_polyN.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       3
#: Number of Node Sets:          10
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
odb = session.odbs['C:/temp/polyN_optimization/execution/UT_EBT_polyN.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Stress components: S11 at Element 1 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xy2 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Stress components: S12 at Element 1 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c2 = session.Curve(xyData=xy2)
xy3 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Stress components: S22 at Element 1 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c3 = session.Curve(xyData=xy3)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, c2, c3, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
session.xyPlots[session.viewports['Viewport: 1'].displayedObject.name].setValues(
    transform=(0.591899, 0, 0, -0.281686, 0, 0.839619, 0, -0.0200252, 0, 0, 
    0.591899, 0, 0, 0, 0, 1))
session.xyPlots['XYPlot-1'].setValues(transform=(8.01348, 0, 0, 0.468854, 0, 
    1.18904, 0, 0.0885544, 0, 0, 0.591899, 0, 0, 0, 0, 1))
session.xyPlots[session.viewports['Viewport: 1'].displayedObject.name].setValues(
    transform=(1.97916, 0, 0, 0.567563, 0, 0.293668, 0, -0.190878, 0, 0, 
    0.146186, 0, 0, 0, 0, 1))
session.xyPlots['XYPlot-1'].setValues(transform=(4.83344, 0, 0, 0.413534, 0, 
    0.614741, 0, 0.579392, 0, 0, 0.146186, 0, 0, 0, 0, 1))
session.xyPlots[session.viewports['Viewport: 1'].displayedObject.name].setValues(
    transform=(0.523155, 0, 0, -0.00692213, 0, 0.0665375, 0, 0.584069, 0, 0, 
    0.0158227, 0, 0, 0, 0, 1))
session.xyPlots['XYPlot-1'].setValues(transform=(0.860524, 0, 0, -0.00183809, 
    0, 0.222059, 0, 0.509729, 0, 0, 0.0158227, 0, 0, 0, 0, 1))
session.xyPlots['XYPlot-1'].setValues(transform=(0.929047, 0, 0, -0.0011257, 0, 
    0.661687, 0, -0.0330147, 0, 0, 0.0158227, 0, 0, 0, 0, 1))
