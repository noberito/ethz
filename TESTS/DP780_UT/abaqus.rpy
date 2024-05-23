# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023 replay file
# Internal Version: 2022_09_28-20.11.55 183150
# Run by cnober on Thu May 23 12:14:27 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=375.819793701172, 
    height=268.650024414062)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='UT_polyN.odb')
#: Model: C:/temp/TESTS/DP780_UT/UT_polyN.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       6
#: Number of Node Sets:          5
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
odb = session.odbs['C:/temp/TESTS/DP780_UT/UT_polyN.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Stress components: S22 at Element 1 Int Point 1 in ELSET ELALL', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
odb = session.odbs['C:/temp/TESTS/DP780_UT/UT_polyN.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Solution dependent state variables: SDV1 at Element 1 Int Point 1 in ELSET ELALL', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xy2 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Stress components: S22 at Element 1 Int Point 1 in ELSET ELALL', 
    steps=('Step-1', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c2 = session.Curve(xyData=xy2)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, c2, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
x0 = session.xyDataObjects['_temp_2']
x1 = session.xyDataObjects['_temp_3']
session.xyReportOptions.setValues(numDigits=8, numberFormat=SCIENTIFIC)
session.writeXYReport(fileName='abaqus.rpt', appendMode=OFF, xyData=(x0, x1))
#* The file designated for the XY Report cannot be opened. Please check that it 
#* is not being used by another process.
x0 = session.xyDataObjects['_temp_2']
x1 = session.xyDataObjects['_temp_3']
session.writeXYReport(fileName='abaqus.rpt', appendMode=OFF, xyData=(x0, x1))
odb = session.odbs['C:/temp/TESTS/DP780_UT/UT_polyN.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
