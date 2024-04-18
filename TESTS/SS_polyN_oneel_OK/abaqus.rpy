# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2023 replay file
# Internal Version: 2022_09_28-20.11.55 183150
# Run by cnober on Tue Apr  9 14:13:16 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=332.922912597656, 
    height=233.010009765625)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='SS_polyN_oneel.odb')
#: Model: C:/temp/TESTS/SS_polyN_oneel/SS_polyN_oneel.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       2
#: Number of Node Sets:          5
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].animationController.setValues(
    animationType=TIME_HISTORY)
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
session.viewports['Viewport: 1'].animationController.stop()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.incrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.decrementFrame()
session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
