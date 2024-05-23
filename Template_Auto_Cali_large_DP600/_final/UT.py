# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.9-1 replay file
# Internal Version: 2009_04_16-15.31.05 92314
# Run by mohr on Thu Aug 04 21:36:02 2011
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=300.00000000, 
    height=200.00000000)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o1 = session.openOdb(name='./temp.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: ./temp.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             ??
#: Number of Element Sets:       ??
#: Number of Node Sets:          ??
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
odb = session.odbs['./temp.odb']
xy0 = session.XYDataFromHistory(name='SDV_EQPS', odb=odb, 
    outputVariableName='Equivalent Plastic Strain: SDV_EQPS at Element 1 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
c0 = session.Curve(xyData=xy0)
xy1 = session.XYDataFromHistory(name='SDV_Seq', odb=odb, 
    outputVariableName='Equivalent stress: SDV_Seq at Element 1 Int Point 1 in ELSET ONSET', 
	steps=('Step-1', ), )
c1 = session.Curve(xyData=xy1)
xy2 = session.XYDataFromHistory(name='SDV_Qeq', odb=odb, 
    outputVariableName='Equivalent Hill stress: SDV_Qeq at Element 1 Int Point 1 in ELSET ONSET', 
	steps=('Step-1', ), )
c2 = session.Curve(xyData=xy1)
session.XYDataFromHistory(name='MISES', odb=odb, 
	outputVariableName='Mises equivalent stress: MISES at Element 1 Int Point 1 in ELSET ONSET', 
	steps=('Step-1', ), )
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c0, c1, ), )
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
x0 = session.xyDataObjects['SDV_EQPS']
x1 = session.xyDataObjects['SDV_Seq']
x2 = session.xyDataObjects['SDV_Qeq']
x3 = session.xyDataObjects['MISES']
session.writeXYReport(fileName='temp.out', appendMode=OFF, xyData=(x0, x1, x2, x3))
