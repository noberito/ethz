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
session.XYDataFromHistory(name='SDV_EQPSdot', odb=odb, 
    outputVariableName='Equivalent Plastic Strain rate: SDV_EQPSdot at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_EQPS', odb=odb, 
    outputVariableName='Equivalent Plastic Strain: SDV_EQPS at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_Seq', odb=odb, 
    outputVariableName='Equivalent stress: SDV_Seq at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_LODE', odb=odb, 
    outputVariableName='Lode parameter: SDV_LODE at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='LEP1', odb=odb, 
    outputVariableName='Principal logarithmic strains: LEP1 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='LEP2', odb=odb, 
    outputVariableName='Principal logarithmic strains: LEP2 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='LEP3', odb=odb, 
    outputVariableName='Principal logarithmic strains: LEP3 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SP1', odb=odb, 
    outputVariableName='Principal stresses: SP1 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SP2', odb=odb, 
    outputVariableName='Principal stresses: SP2 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SP3', odb=odb, 
    outputVariableName='Principal stresses: SP3 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='RF2', odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 4923 in NSET BASE', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='RF2L', odb=odb, 
    outputVariableName='Reaction force: RF2 at Node 1683 in NSET LOAD', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='U2', odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 1683 in NSET LOAD', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='S11', odb=odb, 
    outputVariableName='Stress components: S11 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='S12', odb=odb, 
    outputVariableName='Stress components: S12 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='S13', odb=odb, 
    outputVariableName='Stress components: S13 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='S22', odb=odb, 
    outputVariableName='Stress components: S22 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='S23', odb=odb, 
    outputVariableName='Stress components: S23 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='S33', odb=odb, 
    outputVariableName='Stress components: S33 at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_Svm', odb=odb, 
	outputVariableName='Mises equivalent stress: MISES at Element 5760 Int Point 1 in ELSET ONSET', 
	steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_T', odb=odb, 
    outputVariableName='Temperature: SDV_T at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_TRIAX', odb=odb, 
    outputVariableName='Triaxiality: SDV_TRIAX at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='SDV_D', odb=odb, 
    outputVariableName='Damage: SDV_D at Element 5760 Int Point 1 in ELSET ONSET', 
    steps=('Step-1', ), )
session.XYDataFromHistory(name='U2_EXT_Y05', odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 3689 in NSET EXT_Y05', 
    steps=('Step-1', ), )

x0 = session.xyDataObjects['S11']
x1 = session.xyDataObjects['S12']
x2 = session.xyDataObjects['S13']
x3 = session.xyDataObjects['S22']
x4 = session.xyDataObjects['S23']
x5 = session.xyDataObjects['S33']
x6 = session.xyDataObjects['SP1']
x7 = session.xyDataObjects['SP2']
x8 = session.xyDataObjects['SP3']
x9 = session.xyDataObjects['SDV_Seq']
x10= session.xyDataObjects['SDV_Svm']
x11= session.xyDataObjects['SDV_D']
x12= session.xyDataObjects['SDV_TRIAX']
x13= session.xyDataObjects['SDV_LODE']
x14= session.xyDataObjects['SDV_EQPS']
x15= session.xyDataObjects['SDV_EQPSdot']
x16= session.xyDataObjects['SDV_T']
x17= session.xyDataObjects['U2']
x18= session.xyDataObjects['RF2']
x19= session.xyDataObjects['RF2L']
x20= session.xyDataObjects['U2_EXT_Y05']
#x21= session.xyDataObjects['LEP1']
#x22= session.xyDataObjects['LEP2']
#x23= session.xyDataObjects['LEP3']

session.writeXYReport(fileName='temp.out', appendMode=OFF, xyData=(x17, x18, x20))
# session.writeXYReport(fileName='temp.out', appendMode=OFF, xyData=(x0, x1, 
    # x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, 
    # x18, x19, x20))