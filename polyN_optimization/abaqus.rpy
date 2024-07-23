# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2023 replay file
# Internal Version: 2022_09_28-20.11.55 183150
# Run by cnober on Fri Jul 19 17:54:54 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=297.310424804688, 
    height=142.020004272461)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o2 = session.openOdb(name='NT20_00_UMAT_10_10.odb')
#: Model: C:/temp/polyN_optimization/NT20_00_UMAT_10_10.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       3
#: Number of Node Sets:          10
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o2)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
#: 
#: Node: PART-1-1.9
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00060e+01,  1.50090e+01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.95471e+00,  1.60443e+01,  7.47884e-01,      -      
#: Deformed coordinates (scaled):     9.95471e+00,  1.60443e+01,  7.47884e-01,      -      
#: Displacement (unscaled):          -5.12879e-02,  1.03525e+00,  8.38482e-05,  1.03652e+00
#: 
#: Node: PART-1-1.11
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00060e+01,  1.50090e+01,  0.00000e+00,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.95479e+00,  1.60443e+01,  0.00000e+00,      -      
#: Deformed coordinates (scaled):     9.95479e+00,  1.60443e+01,  0.00000e+00,      -      
#: Displacement (unscaled):          -5.12084e-02,  1.03525e+00, -0.00000e+00,  1.03652e+00
#: 
#: Nodes for distance: PART-1-1.9, PART-1-1.11
#:                                        1             2             3        Magnitude
#: Base distance:                     0.00000e+00,  0.00000e+00, -7.47800e-01,  7.47800e-01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):      7.91550e-05,  0.00000e+00, -7.47884e-01,  7.47884e-01
#: Deformed distance (scaled):        7.91550e-05,  0.00000e+00, -7.47884e-01,  7.47884e-01
#: Relative displacement (unscaled):  7.95312e-05,  0.00000e+00, -8.38482e-05,  1.15567e-04
#: 
#: Node: PART-1-1.10
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  1.50090e+01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   0.00000e+00,  1.60443e+01,  7.45683e-01,      -      
#: Deformed coordinates (scaled):     0.00000e+00,  1.60443e+01,  7.45683e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  1.03525e+00, -2.11707e-03,  1.03525e+00
#: 
#: Node: PART-1-1.9
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.00060e+01,  1.50090e+01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   9.95471e+00,  1.60443e+01,  7.47884e-01,      -      
#: Deformed coordinates (scaled):     9.95471e+00,  1.60443e+01,  7.47884e-01,      -      
#: Displacement (unscaled):          -5.12879e-02,  1.03525e+00,  8.38482e-05,  1.03652e+00
#: 
#: Nodes for distance: PART-1-1.10, PART-1-1.9
#:                                        1             2             3        Magnitude
#: Base distance:                     1.00060e+01,  0.00000e+00,  0.00000e+00,  1.00060e+01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):      9.95471e+00,  0.00000e+00,  2.20090e-03,  9.95471e+00
#: Deformed distance (scaled):        9.95471e+00,  0.00000e+00,  2.20090e-03,  9.95471e+00
#: Relative displacement (unscaled): -5.12879e-02,  0.00000e+00,  2.20091e-03,  5.13351e-02
#: 
#: Node: PART-1-1.10
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  1.50090e+01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   0.00000e+00,  1.60443e+01,  7.45683e-01,      -      
#: Deformed coordinates (scaled):     0.00000e+00,  1.60443e+01,  7.45683e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  1.03525e+00, -2.11707e-03,  1.03525e+00
#: 
#: Node: PART-1-1.4136
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.25178e-01,  4.00243e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.15868e-01,  5.84297e-01,  5.59323e-01,      -      
#: Deformed coordinates (scaled):     1.15868e-01,  5.84297e-01,  5.59323e-01,      -      
#: Displacement (unscaled):          -9.31012e-03,  1.84054e-01, -1.88477e-01,  2.63602e-01
#: 
#: Nodes for distance: PART-1-1.10, PART-1-1.4136
#:                                        1             2             3        Magnitude
#: Base distance:                     1.25178e-01, -1.46088e+01,  0.00000e+00,  1.46093e+01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):      1.15868e-01, -1.54600e+01, -1.86360e-01,  1.54615e+01
#: Deformed distance (scaled):        1.15868e-01, -1.54600e+01, -1.86360e-01,  1.54615e+01
#: Relative displacement (unscaled): -9.31012e-03, -8.51196e-01, -1.86360e-01,  8.71408e-01
#: 
#: Node: PART-1-1.4134
#:                                         1             2             3        Magnitude
#: Base coordinates:                  1.25103e-01,  2.00122e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   1.15250e-01,  2.95070e-01,  5.37489e-01,      -      
#: Deformed coordinates (scaled):     1.15250e-01,  2.95070e-01,  5.37489e-01,      -      
#: Displacement (unscaled):          -9.85262e-03,  9.49477e-02, -2.10311e-01,  2.30961e-01
#: 
#: Node: PART-1-1.3685
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  1.00060e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   7.43600e-37,  1.47978e-01,  5.30756e-01,      -      
#: Deformed coordinates (scaled):     7.43600e-37,  1.47978e-01,  5.30756e-01,      -      
#: Displacement (unscaled):           7.43600e-37,  4.79178e-02, -2.17043e-01,  2.22270e-01
#: 
#: Nodes for distance: PART-1-1.4134, PART-1-1.3685
#:                                        1             2             3        Magnitude
#: Base distance:                    -1.25103e-01, -1.00062e-01,  0.00000e+00,  1.60197e-01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):     -1.15250e-01, -1.47092e-01, -6.73205e-03,  1.86987e-01
#: Deformed distance (scaled):       -1.15250e-01, -1.47092e-01, -6.73205e-03,  1.86987e-01
#: Relative displacement (unscaled):  9.85262e-03, -4.70299e-02, -6.73200e-03,  4.85201e-02
#: 
#: Node: PART-1-1.3685
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  1.00060e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   7.43600e-37,  1.47978e-01,  5.30756e-01,      -      
#: Deformed coordinates (scaled):     7.43600e-37,  1.47978e-01,  5.30756e-01,      -      
#: Displacement (unscaled):           7.43600e-37,  4.79178e-02, -2.17043e-01,  2.22270e-01
#: 
#: Node: PART-1-1.10
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  1.50090e+01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   0.00000e+00,  1.60443e+01,  7.45683e-01,      -      
#: Deformed coordinates (scaled):     0.00000e+00,  1.60443e+01,  7.45683e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  1.03525e+00, -2.11707e-03,  1.03525e+00
#: 
#: Nodes for distance: PART-1-1.3685, PART-1-1.10
#:                                        1             2             3        Magnitude
#: Base distance:                     0.00000e+00,  1.49089e+01,  0.00000e+00,  1.49089e+01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):     -7.43600e-37,  1.58963e+01,  2.14926e-01,  1.58977e+01
#: Deformed distance (scaled):       -7.43600e-37,  1.58963e+01,  2.14926e-01,  1.58977e+01
#: Relative displacement (unscaled): -7.43600e-37,  9.87332e-01,  2.14926e-01,  1.01045e+00
session.viewports['Viewport: 1'].view.setValues(nearPlane=24.5461, 
    farPlane=48.9614, width=34.932, height=15.7382, viewOffsetX=5.16632, 
    viewOffsetY=-1.51229)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=0 )
#: 
#: Node: PART-1-1.8675
#:                                         1             2             3        Magnitude
#: Base coordinates:                  4.64307e+00,  5.04048e+00,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   4.64307e+00,  5.04048e+00,  7.47800e-01,      -      
#: Deformed coordinates (scaled):     4.64307e+00,  5.04048e+00,  7.47800e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: PART-1-1.3873
#:                                         1             2             3        Magnitude
#: Base coordinates:                  2.50606e-01,  6.00368e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   2.50606e-01,  6.00368e-01,  7.47800e-01,      -      
#: Deformed coordinates (scaled):     2.50606e-01,  6.00368e-01,  7.47800e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Nodes for distance: PART-1-1.8675, PART-1-1.3873
#:                                        1             2             3        Magnitude
#: Base distance:                    -4.39246e+00, -4.44012e+00,  0.00000e+00,  6.24567e+00
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):     -4.39246e+00, -4.44012e+00,  0.00000e+00,  6.24567e+00
#: Deformed distance (scaled):       -4.39246e+00, -4.44012e+00,  0.00000e+00,  6.24567e+00
#: Relative displacement (unscaled):  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: PART-1-1.3869
#:                                         1             2             3        Magnitude
#: Base coordinates:                  2.50206e-01,  2.00124e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   2.50206e-01,  2.00124e-01,  7.47800e-01,      -      
#: Deformed coordinates (scaled):     2.50206e-01,  2.00124e-01,  7.47800e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: PART-1-1.253
#:                                         1             2             3        Magnitude
#: Base coordinates:                  2.50150e-01,  1.50090e+01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   2.50150e-01,  1.50090e+01,  7.47800e-01,      -      
#: Deformed coordinates (scaled):     2.50150e-01,  1.50090e+01,  7.47800e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Nodes for distance: PART-1-1.3869, PART-1-1.253
#:                                        1             2             3        Magnitude
#: Base distance:                    -5.59390e-05,  1.48089e+01,  0.00000e+00,  1.48089e+01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):     -5.59390e-05,  1.48089e+01,  0.00000e+00,  1.48089e+01
#: Deformed distance (scaled):       -5.59390e-05,  1.48089e+01,  0.00000e+00,  1.48089e+01
#: Relative displacement (unscaled):  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=33.9895, 
    farPlane=38.2268, width=18.7791, height=10.3411, viewOffsetX=-1.73346, 
    viewOffsetY=-1.97851)
session.viewports['Viewport: 1'].view.setValues(nearPlane=34.0843, 
    farPlane=38.132, width=18.8315, height=10.37, viewOffsetX=-3.16351, 
    viewOffsetY=-7.44376)
#: 
#: Node: PART-1-1.3689
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  5.00300e-01,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   0.00000e+00,  5.00300e-01,  7.47800e-01,      -      
#: Deformed coordinates (scaled):     0.00000e+00,  5.00300e-01,  7.47800e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Node: PART-1-1.3448
#:                                         1             2             3        Magnitude
#: Base coordinates:                  0.00000e+00,  0.00000e+00,  7.47800e-01,      -      
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed coordinates (unscaled):   0.00000e+00,  0.00000e+00,  7.47800e-01,      -      
#: Deformed coordinates (scaled):     0.00000e+00,  0.00000e+00,  7.47800e-01,      -      
#: Displacement (unscaled):           0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
#: 
#: Nodes for distance: PART-1-1.3689, PART-1-1.3448
#:                                        1             2             3        Magnitude
#: Base distance:                     0.00000e+00, -5.00300e-01,  0.00000e+00,  5.00300e-01
#: Scale:                             1.00000e+00,  1.00000e+00,  1.00000e+00,      -      
#: Deformed distance (unscaled):      0.00000e+00, -5.00300e-01,  0.00000e+00,  5.00300e-01
#: Deformed distance (scaled):        0.00000e+00, -5.00300e-01,  0.00000e+00,  5.00300e-01
#: Relative displacement (unscaled):  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00
