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
xy0 = session.XYDataFromHistory(name='U2', odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 25791 in NSET LOAD', 
    steps=('Step-1', ), )
xy1 = session.XYDataFromHistory(name='RF2', odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 14364 in NSET BASE', steps=(
    'Step-1', ), )
xy2 =session.XYDataFromHistory(name='EXT_X05_Y05', odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 13185 in NSET EXT_X05_Y05', 
    steps=('Step-1', ), )	
x0 = session.xyDataObjects['U2']
x1 = session.xyDataObjects['RF2']
x2 = session.xyDataObjects['EXT_X05_Y05']
session.writeXYReport(fileName='temp.out', appendMode=OFF, xyData=(x0, x1, x2))