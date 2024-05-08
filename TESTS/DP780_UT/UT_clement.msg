1

   Abaqus 2023                                  Date 07-mai-2024   Time 14:11:53
   For use by ETH ZUERICH under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                2.000E-02
          AND A TOTAL TIME PERIOD OF                            1.00    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-05
          THE MAXIMUM TIME INCREMENT ALLOWED IS                3.000E-02

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE FORCE                               1.000E-02
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                  5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     100

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH      5.00    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.000E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       599.       TIME AVG. FORCE        599.    
 LARGEST RESIDUAL FORCE             -231.       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       515.       TIME AVG. FORCE        515.    
 LARGEST RESIDUAL FORCE             -90.8       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.904E-03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       459.       TIME AVG. FORCE        459.    
 LARGEST RESIDUAL FORCE             -35.4       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.168E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       441.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE             -14.6       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.138E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       434.       TIME AVG. FORCE        434.    
 LARGEST RESIDUAL FORCE             -6.03       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.692E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       430.       TIME AVG. FORCE        430.    
 LARGEST RESIDUAL FORCE             -2.51       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.596E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       429.       TIME AVG. FORCE        429.    
 LARGEST RESIDUAL FORCE             -1.05       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.000E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.495E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.000E-02,  FRACTION OF STEP COMPLETED  2.000E-02
 STEP TIME COMPLETED       2.000E-02,  TOTAL TIME COMPLETED        2.000E-02


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.000E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.038E+03   TIME AVG. FORCE       2.733E+03
 LARGEST RESIDUAL FORCE             5.949E+03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -8.531E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.504E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.088E+05   TIME AVG. FORCE       1.546E+05
 LARGEST RESIDUAL FORCE            -3.338E+05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          2.90       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         2.99       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.037E+06   TIME AVG. FORCE       1.519E+06
 LARGEST RESIDUAL FORCE            -2.519E+06   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          10.8       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         7.94       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0
 

 ***WARNING: THE SYSTEM MATRIX HAS 3 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.422E+07   TIME AVG. FORCE       7.112E+06
 LARGEST RESIDUAL FORCE            -8.078E+06   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          26.7       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         15.9       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     2 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  5.000E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       486.       TIME AVG. FORCE        458.    
 LARGEST RESIDUAL FORCE             -59.5       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.000E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.994E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       450.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE             -21.9       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.000E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.147E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       438.       TIME AVG. FORCE        434.    
 LARGEST RESIDUAL FORCE             -8.43       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.015E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.386E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       433.       TIME AVG. FORCE        431.    
 LARGEST RESIDUAL FORCE             -3.24       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.068E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.288E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       432.       TIME AVG. FORCE        430.    
 LARGEST RESIDUAL FORCE             -1.25       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.089E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.035E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.000E-03,  FRACTION OF STEP COMPLETED  2.500E-02
 STEP TIME COMPLETED       2.500E-02,  TOTAL TIME COMPLETED        2.500E-02


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.000E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       433.       TIME AVG. FORCE        431.    
 LARGEST RESIDUAL FORCE            -0.166       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.092E-03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.897E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.000E-03,  FRACTION OF STEP COMPLETED  3.000E-02
 STEP TIME COMPLETED       3.000E-02,  TOTAL TIME COMPLETED        3.000E-02


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.000E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       435.       TIME AVG. FORCE        432.    
 LARGEST RESIDUAL FORCE             0.733       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.079E-03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.246E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.500E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.000E-03,  FRACTION OF STEP COMPLETED  3.500E-02
 STEP TIME COMPLETED       3.500E-02,  TOTAL TIME COMPLETED        3.500E-02


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.500E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       439.       TIME AVG. FORCE        434.    
 LARGEST RESIDUAL FORCE              2.27       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -7.602E-03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.699E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       439.       TIME AVG. FORCE        434.    
 LARGEST RESIDUAL FORCE              1.30       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -7.592E-03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.680E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.125E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.500E-03,  FRACTION OF STEP COMPLETED  4.250E-02
 STEP TIME COMPLETED       4.250E-02,  TOTAL TIME COMPLETED        4.250E-02


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       445.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              5.69       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.138E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.028E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              4.92       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.137E-02   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.086E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              4.31       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.137E-02   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.159E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              3.76       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.136E-02   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.445E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              3.28       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.136E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.812E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              2.87       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.135E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.172E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              2.50       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.135E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.638E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       444.       TIME AVG. FORCE        435.    
 LARGEST RESIDUAL FORCE              2.16       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.135E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.151E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  5.375E-02
 STEP TIME COMPLETED       5.375E-02,  TOTAL TIME COMPLETED        5.375E-02


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       450.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              7.08       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.133E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.243E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              6.06       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.132E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.883E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              5.22       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.132E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.438E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              4.44       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.131E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.317E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              3.84       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.130E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.242E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              3.28       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.130E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.418E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              2.79       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.129E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.562E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              2.38       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.129E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.895E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       449.       TIME AVG. FORCE        437.    
 LARGEST RESIDUAL FORCE              2.06       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.129E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.322E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   9 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  9 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  6.500E-02
 STEP TIME COMPLETED       6.500E-02,  TOTAL TIME COMPLETED        6.500E-02


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              5.95       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.128E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.018E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              5.00       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.127E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.098E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              4.24       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.126E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.609E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              3.57       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.125E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.509E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              3.00       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.536E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              2.56       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.587E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       454.       TIME AVG. FORCE        439.    
 LARGEST RESIDUAL FORCE              2.16       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.906E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  7.625E-02
 STEP TIME COMPLETED       7.625E-02,  TOTAL TIME COMPLETED        7.625E-02


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       459.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              6.11       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.159E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       459.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              5.05       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.000E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       459.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              4.22       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.304E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       458.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              3.49       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.901E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       458.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              2.91       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.718E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       458.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              2.40       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.790E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       458.       TIME AVG. FORCE        441.    
 LARGEST RESIDUAL FORCE              2.00       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.989E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  8.750E-02
 STEP TIME COMPLETED       8.750E-02,  TOTAL TIME COMPLETED        8.750E-02


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       464.       TIME AVG. FORCE        444.    
 LARGEST RESIDUAL FORCE              7.73       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.148E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       463.       TIME AVG. FORCE        444.    
 LARGEST RESIDUAL FORCE              4.47       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.293E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       463.       TIME AVG. FORCE        444.    
 LARGEST RESIDUAL FORCE              3.54       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.288E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       463.       TIME AVG. FORCE        444.    
 LARGEST RESIDUAL FORCE              2.89       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.953E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       463.       TIME AVG. FORCE        444.    
 LARGEST RESIDUAL FORCE              2.36       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.274E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       463.       TIME AVG. FORCE        444.    
 LARGEST RESIDUAL FORCE              1.93       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.301E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  9.875E-02
 STEP TIME COMPLETED       9.875E-02,  TOTAL TIME COMPLETED        9.875E-02


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       468.       TIME AVG. FORCE        446.    
 LARGEST RESIDUAL FORCE              5.56       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.226E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       468.       TIME AVG. FORCE        446.    
 LARGEST RESIDUAL FORCE              4.47       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.034E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       467.       TIME AVG. FORCE        446.    
 LARGEST RESIDUAL FORCE              3.61       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.391E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       467.       TIME AVG. FORCE        446.    
 LARGEST RESIDUAL FORCE              2.92       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.755E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       467.       TIME AVG. FORCE        446.    
 LARGEST RESIDUAL FORCE              2.36       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.474E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       467.       TIME AVG. FORCE        446.    
 LARGEST RESIDUAL FORCE              1.89       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.506E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.110    
 STEP TIME COMPLETED       0.110    ,  TOTAL TIME COMPLETED        0.110    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       472.       TIME AVG. FORCE        448.    
 LARGEST RESIDUAL FORCE              5.46       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.325E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       472.       TIME AVG. FORCE        448.    
 LARGEST RESIDUAL FORCE              4.37       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.074E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       472.       TIME AVG. FORCE        448.    
 LARGEST RESIDUAL FORCE              3.46       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.626E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       472.       TIME AVG. FORCE        448.    
 LARGEST RESIDUAL FORCE              2.75       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.818E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       471.       TIME AVG. FORCE        448.    
 LARGEST RESIDUAL FORCE              2.17       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.442E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.121    
 STEP TIME COMPLETED       0.121    ,  TOTAL TIME COMPLETED        0.121    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       476.       TIME AVG. FORCE        450.    
 LARGEST RESIDUAL FORCE              5.91       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.533E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       476.       TIME AVG. FORCE        450.    
 LARGEST RESIDUAL FORCE              4.60       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.227E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       476.       TIME AVG. FORCE        450.    
 LARGEST RESIDUAL FORCE              3.63       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.515E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       476.       TIME AVG. FORCE        450.    
 LARGEST RESIDUAL FORCE              2.85       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.505E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       476.       TIME AVG. FORCE        450.    
 LARGEST RESIDUAL FORCE              2.21       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.925E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.132    
 STEP TIME COMPLETED       0.132    ,  TOTAL TIME COMPLETED        0.132    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       480.       TIME AVG. FORCE        452.    
 LARGEST RESIDUAL FORCE              5.61       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.563E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       480.       TIME AVG. FORCE        452.    
 LARGEST RESIDUAL FORCE              4.34       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.212E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       480.       TIME AVG. FORCE        452.    
 LARGEST RESIDUAL FORCE              3.36       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.373E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       480.       TIME AVG. FORCE        452.    
 LARGEST RESIDUAL FORCE              2.60       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.261E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       480.       TIME AVG. FORCE        452.    
 LARGEST RESIDUAL FORCE              2.02       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.606E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.144    
 STEP TIME COMPLETED       0.144    ,  TOTAL TIME COMPLETED        0.144    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       484.       TIME AVG. FORCE        454.    
 LARGEST RESIDUAL FORCE              5.20       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.509E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       484.       TIME AVG. FORCE        454.    
 LARGEST RESIDUAL FORCE              3.97       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.172E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       484.       TIME AVG. FORCE        454.    
 LARGEST RESIDUAL FORCE              3.02       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.938E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       484.       TIME AVG. FORCE        454.    
 LARGEST RESIDUAL FORCE              2.30       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.802E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       484.       TIME AVG. FORCE        454.    
 LARGEST RESIDUAL FORCE              1.75       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.173E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.155    
 STEP TIME COMPLETED       0.155    ,  TOTAL TIME COMPLETED        0.155    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       488.       TIME AVG. FORCE        456.    
 LARGEST RESIDUAL FORCE              4.83       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.487E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       488.       TIME AVG. FORCE        456.    
 LARGEST RESIDUAL FORCE              3.64       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.129E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       488.       TIME AVG. FORCE        456.    
 LARGEST RESIDUAL FORCE              2.76       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.517E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       488.       TIME AVG. FORCE        456.    
 LARGEST RESIDUAL FORCE              2.08       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.484E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.166    
 STEP TIME COMPLETED       0.166    ,  TOTAL TIME COMPLETED        0.166    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.125E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       492.       TIME AVG. FORCE        458.    
 LARGEST RESIDUAL FORCE              5.52       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.843E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       492.       TIME AVG. FORCE        458.    
 LARGEST RESIDUAL FORCE              4.09       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.340E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       492.       TIME AVG. FORCE        458.    
 LARGEST RESIDUAL FORCE              3.05       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.878E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       491.       TIME AVG. FORCE        458.    
 LARGEST RESIDUAL FORCE              2.25       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.125E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.390E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.688E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.125E-02,  FRACTION OF STEP COMPLETED  0.178    
 STEP TIME COMPLETED       0.178    ,  TOTAL TIME COMPLETED        0.178    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.688E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       535.       TIME AVG. FORCE        463.    
 LARGEST RESIDUAL FORCE              238.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.212E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       534.       TIME AVG. FORCE        463.    
 LARGEST RESIDUAL FORCE              234.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.732E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       533.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              228.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.102E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       532.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              221.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.572E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       530.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              212.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.170E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       529.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              203.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.829E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       527.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              193.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.533E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       526.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              183.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.425E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       524.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              173.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.688E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.330E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    18 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  8.438E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       494.       TIME AVG. FORCE        460.    
 LARGEST RESIDUAL FORCE              3.04       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.438E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.506E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       494.       TIME AVG. FORCE        460.    
 LARGEST RESIDUAL FORCE              1.65       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.438E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.322E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.266E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.438E-03,  FRACTION OF STEP COMPLETED  0.186    
 STEP TIME COMPLETED       0.186    ,  TOTAL TIME COMPLETED        0.186    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.266E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       499.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              4.89       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.410E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       499.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              4.03       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.000E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       499.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              3.32       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.514E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       498.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              2.74       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.369E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       498.       TIME AVG. FORCE        462.    
 LARGEST RESIDUAL FORCE              2.27       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.530E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.266E-02,  FRACTION OF STEP COMPLETED  0.199    
 STEP TIME COMPLETED       0.199    ,  TOTAL TIME COMPLETED        0.199    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.266E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       503.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              7.56       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.939E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       503.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              6.14       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.308E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       503.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              4.92       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.055E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       503.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              4.10       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.683E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       503.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              3.32       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.307E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       502.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              2.70       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.702E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       502.       TIME AVG. FORCE        464.    
 LARGEST RESIDUAL FORCE              2.22       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.656E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.266E-02,  FRACTION OF STEP COMPLETED  0.211    
 STEP TIME COMPLETED       0.211    ,  TOTAL TIME COMPLETED        0.211    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.266E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       507.       TIME AVG. FORCE        466.    
 LARGEST RESIDUAL FORCE              6.62       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.177E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       507.       TIME AVG. FORCE        466.    
 LARGEST RESIDUAL FORCE              5.22       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.422E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       507.       TIME AVG. FORCE        466.    
 LARGEST RESIDUAL FORCE              4.17       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.980E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       506.       TIME AVG. FORCE        466.    
 LARGEST RESIDUAL FORCE              3.35       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.683E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       506.       TIME AVG. FORCE        466.    
 LARGEST RESIDUAL FORCE              2.68       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.100E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       506.       TIME AVG. FORCE        466.    
 LARGEST RESIDUAL FORCE              2.16       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.865E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.266E-02,  FRACTION OF STEP COMPLETED  0.224    
 STEP TIME COMPLETED       0.224    ,  TOTAL TIME COMPLETED        0.224    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.266E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       511.       TIME AVG. FORCE        468.    
 LARGEST RESIDUAL FORCE              7.06       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.107E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       511.       TIME AVG. FORCE        468.    
 LARGEST RESIDUAL FORCE              5.21       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.390E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       510.       TIME AVG. FORCE        468.    
 LARGEST RESIDUAL FORCE              4.13       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.035E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       510.       TIME AVG. FORCE        468.    
 LARGEST RESIDUAL FORCE              3.27       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.024E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       510.       TIME AVG. FORCE        468.    
 LARGEST RESIDUAL FORCE              2.59       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.243E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       510.       TIME AVG. FORCE        468.    
 LARGEST RESIDUAL FORCE              2.05       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.963E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.266E-02,  FRACTION OF STEP COMPLETED  0.237    
 STEP TIME COMPLETED       0.237    ,  TOTAL TIME COMPLETED        0.237    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.266E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       515.       TIME AVG. FORCE        470.    
 LARGEST RESIDUAL FORCE              8.53       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.377E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       514.       TIME AVG. FORCE        470.    
 LARGEST RESIDUAL FORCE              4.87       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.131E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       514.       TIME AVG. FORCE        470.    
 LARGEST RESIDUAL FORCE              3.69       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.190E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       514.       TIME AVG. FORCE        470.    
 LARGEST RESIDUAL FORCE              2.89       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.664E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       514.       TIME AVG. FORCE        470.    
 LARGEST RESIDUAL FORCE              4.24       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.266E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.776E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    23 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.164E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       511.       TIME AVG. FORCE        470.    
 LARGEST RESIDUAL FORCE             0.558       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.164E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.451E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.164E-03,  FRACTION OF STEP COMPLETED  0.240    
 STEP TIME COMPLETED       0.240    ,  TOTAL TIME COMPLETED        0.240    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.164E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       512.       TIME AVG. FORCE        472.    
 LARGEST RESIDUAL FORCE            -0.108       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.164E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.975E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.746E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.164E-03,  FRACTION OF STEP COMPLETED  0.243    
 STEP TIME COMPLETED       0.243    ,  TOTAL TIME COMPLETED        0.243    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.746E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       514.       TIME AVG. FORCE        474.    
 LARGEST RESIDUAL FORCE            -0.124       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.746E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.309E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.119E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.746E-03,  FRACTION OF STEP COMPLETED  0.248    
 STEP TIME COMPLETED       0.248    ,  TOTAL TIME COMPLETED        0.248    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.119E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       516.       TIME AVG. FORCE        475.    
 LARGEST RESIDUAL FORCE             0.633       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.119E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.038E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.068E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.119E-03,  FRACTION OF STEP COMPLETED  0.255    
 STEP TIME COMPLETED       0.255    ,  TOTAL TIME COMPLETED        0.255    


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.068E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       519.       TIME AVG. FORCE        477.    
 LARGEST RESIDUAL FORCE              3.57       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.068E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.637E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       519.       TIME AVG. FORCE        477.    
 LARGEST RESIDUAL FORCE              2.30       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.068E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.167E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.602E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.068E-02,  FRACTION OF STEP COMPLETED  0.265    
 STEP TIME COMPLETED       0.265    ,  TOTAL TIME COMPLETED        0.265    


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.602E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       539.       TIME AVG. FORCE        479.    
 LARGEST RESIDUAL FORCE              95.3       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.602E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.559E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       541.       TIME AVG. FORCE        479.    
 LARGEST RESIDUAL FORCE              105.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.602E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.836E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       540.       TIME AVG. FORCE        479.    
 LARGEST RESIDUAL FORCE              103.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.602E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.152E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       530.       TIME AVG. FORCE        479.    
 LARGEST RESIDUAL FORCE              90.5       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.602E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.483E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       543.       TIME AVG. FORCE        479.    
 LARGEST RESIDUAL FORCE              87.3       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.602E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.245E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       539.       TIME AVG. FORCE        479.    
 LARGEST RESIDUAL FORCE              97.7       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.602E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.023E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    28 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  4.005E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       521.       TIME AVG. FORCE        478.    
 LARGEST RESIDUAL FORCE             0.944       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.005E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.152E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    6.007E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.005E-03,  FRACTION OF STEP COMPLETED  0.269    
 STEP TIME COMPLETED       0.269    ,  TOTAL TIME COMPLETED        0.269    


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.007E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       738.       TIME AVG. FORCE        487.    
 LARGEST RESIDUAL FORCE            -1.468E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.007E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.372E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    9.010E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.007E-03,  FRACTION OF STEP COMPLETED  0.275    
 STEP TIME COMPLETED       0.275    ,  TOTAL TIME COMPLETED        0.275    


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.010E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       526.       TIME AVG. FORCE        489.    
 LARGEST RESIDUAL FORCE             0.574       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.010E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.353E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.352E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.010E-03,  FRACTION OF STEP COMPLETED  0.284    
 STEP TIME COMPLETED       0.284    ,  TOTAL TIME COMPLETED        0.284    


  INCREMENT    31 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.352E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       530.       TIME AVG. FORCE        490.    
 LARGEST RESIDUAL FORCE              6.39       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.045E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       530.       TIME AVG. FORCE        490.    
 LARGEST RESIDUAL FORCE              10.5       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.359E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       530.       TIME AVG. FORCE        490.    
 LARGEST RESIDUAL FORCE              3.94       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.136E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       530.       TIME AVG. FORCE        490.    
 LARGEST RESIDUAL FORCE              3.62       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.079E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       531.       TIME AVG. FORCE        490.    
 LARGEST RESIDUAL FORCE              1.20       AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.848E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.352E-02,  FRACTION OF STEP COMPLETED  0.298    
 STEP TIME COMPLETED       0.298    ,  TOTAL TIME COMPLETED        0.298    


  INCREMENT    32 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.352E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       536.       TIME AVG. FORCE        491.    
 LARGEST RESIDUAL FORCE              7.70       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.927E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       536.       TIME AVG. FORCE        491.    
 LARGEST RESIDUAL FORCE              2.61       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.959E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       536.       TIME AVG. FORCE        491.    
 LARGEST RESIDUAL FORCE              6.67       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.227E-06   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       536.       TIME AVG. FORCE        491.    
 LARGEST RESIDUAL FORCE              4.56       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.295E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       536.       TIME AVG. FORCE        491.    
 LARGEST RESIDUAL FORCE              6.21       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.784E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       537.       TIME AVG. FORCE        492.    
 LARGEST RESIDUAL FORCE              1.10       AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.825E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.352E-02,  FRACTION OF STEP COMPLETED  0.311    
 STEP TIME COMPLETED       0.311    ,  TOTAL TIME COMPLETED        0.311    


  INCREMENT    33 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.352E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       543.       TIME AVG. FORCE        493.    
 LARGEST RESIDUAL FORCE              6.15       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.998E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       765.       TIME AVG. FORCE        500.    
 LARGEST RESIDUAL FORCE             -362.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.131E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       564.       TIME AVG. FORCE        494.    
 LARGEST RESIDUAL FORCE             -214.       AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.455E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.331E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       547.       TIME AVG. FORCE        493.    
 LARGEST RESIDUAL FORCE              72.0       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.503E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       544.       TIME AVG. FORCE        493.    
 LARGEST RESIDUAL FORCE              44.3       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.546E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       544.       TIME AVG. FORCE        493.    
 LARGEST RESIDUAL FORCE              51.0       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.122E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       545.       TIME AVG. FORCE        493.    
 LARGEST RESIDUAL FORCE              57.2       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.352E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.424E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    33 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.379E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       538.       TIME AVG. FORCE        493.    
 LARGEST RESIDUAL FORCE             0.405       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.379E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.081E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.379E-03,  FRACTION OF STEP COMPLETED  0.315    
 STEP TIME COMPLETED       0.315    ,  TOTAL TIME COMPLETED        0.315    


  INCREMENT    34 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.379E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       540.       TIME AVG. FORCE        494.    
 LARGEST RESIDUAL FORCE             0.297       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.379E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.427E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.068E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.379E-03,  FRACTION OF STEP COMPLETED  0.318    
 STEP TIME COMPLETED       0.318    ,  TOTAL TIME COMPLETED        0.318    


  INCREMENT    35 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.068E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       542.       TIME AVG. FORCE        496.    
 LARGEST RESIDUAL FORCE             0.246       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.068E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.333E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.602E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.068E-03,  FRACTION OF STEP COMPLETED  0.323    
 STEP TIME COMPLETED       0.323    ,  TOTAL TIME COMPLETED        0.323    


  INCREMENT    36 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.602E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       546.       TIME AVG. FORCE        497.    
 LARGEST RESIDUAL FORCE             0.916       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.602E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.942E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.140E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.602E-03,  FRACTION OF STEP COMPLETED  0.331    
 STEP TIME COMPLETED       0.331    ,  TOTAL TIME COMPLETED        0.331    


  INCREMENT    37 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.140E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       552.       TIME AVG. FORCE        499.    
 LARGEST RESIDUAL FORCE              8.28       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.140E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.765E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       552.       TIME AVG. FORCE        499.    
 LARGEST RESIDUAL FORCE              5.81       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.140E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.565E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       552.       TIME AVG. FORCE        499.    
 LARGEST RESIDUAL FORCE              3.76       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.140E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.142E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       775.       TIME AVG. FORCE        505.    
 LARGEST RESIDUAL FORCE             -334.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.140E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.123E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    37 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  2.851E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       547.       TIME AVG. FORCE        498.    
 LARGEST RESIDUAL FORCE             0.296       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.851E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.667E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.276E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.851E-03,  FRACTION OF STEP COMPLETED  0.334    
 STEP TIME COMPLETED       0.334    ,  TOTAL TIME COMPLETED        0.334    


  INCREMENT    38 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.276E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       550.       TIME AVG. FORCE        500.    
 LARGEST RESIDUAL FORCE             8.045E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.276E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.995E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    6.414E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.276E-03,  FRACTION OF STEP COMPLETED  0.338    
 STEP TIME COMPLETED       0.338    ,  TOTAL TIME COMPLETED        0.338    


  INCREMENT    39 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.414E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       554.       TIME AVG. FORCE        501.    
 LARGEST RESIDUAL FORCE             0.368       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.414E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.346E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    9.622E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.414E-03,  FRACTION OF STEP COMPLETED  0.344    
 STEP TIME COMPLETED       0.344    ,  TOTAL TIME COMPLETED        0.344    


  INCREMENT    40 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.622E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       590.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE              311.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.297E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       629.       TIME AVG. FORCE        504.    
 LARGEST RESIDUAL FORCE             -165.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.438E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       602.       TIME AVG. FORCE        504.    
 LARGEST RESIDUAL FORCE             -75.1       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.104E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       601.       TIME AVG. FORCE        504.    
 LARGEST RESIDUAL FORCE             -74.7       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.116E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       590.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE             -56.1       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.662E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       581.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE             -40.6       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.088E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       581.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE             -41.3       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.797E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       572.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE             -23.6       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.961E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       574.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE             -27.7       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.622E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.695E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    40 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  4.811E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       557.       TIME AVG. FORCE        503.    
 LARGEST RESIDUAL FORCE             0.579       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.811E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.262E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.216E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.811E-03,  FRACTION OF STEP COMPLETED  0.349    
 STEP TIME COMPLETED       0.349    ,  TOTAL TIME COMPLETED        0.349    


  INCREMENT    41 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.216E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       562.       TIME AVG. FORCE        504.    
 LARGEST RESIDUAL FORCE              1.56       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.216E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.968E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.082E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.216E-03,  FRACTION OF STEP COMPLETED  0.356    
 STEP TIME COMPLETED       0.356    ,  TOTAL TIME COMPLETED        0.356    


  INCREMENT    42 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.082E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       571.       TIME AVG. FORCE        506.    
 LARGEST RESIDUAL FORCE              5.21       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.082E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.093E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       571.       TIME AVG. FORCE        506.    
 LARGEST RESIDUAL FORCE              3.51       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.082E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.641E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       571.       TIME AVG. FORCE        506.    
 LARGEST RESIDUAL FORCE              2.35       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.082E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.052E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.624E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.082E-02,  FRACTION OF STEP COMPLETED  0.367    
 STEP TIME COMPLETED       0.367    ,  TOTAL TIME COMPLETED        0.367    


  INCREMENT    43 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.624E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       587.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE              15.1       AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.971E-04   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       586.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE              34.3       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.145E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       594.       TIME AVG. FORCE        508.    
 LARGEST RESIDUAL FORCE              192.       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.571E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       599.       TIME AVG. FORCE        508.    
 LARGEST RESIDUAL FORCE              78.1       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.161E-03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       589.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE              21.7       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.674E-03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       589.       TIME AVG. FORCE        508.    
 LARGEST RESIDUAL FORCE             -8.26       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.867E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       594.       TIME AVG. FORCE        508.    
 LARGEST RESIDUAL FORCE             -15.7       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.004E-04   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       587.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE             -3.12       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.560E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       592.       TIME AVG. FORCE        508.    
 LARGEST RESIDUAL FORCE             -11.3       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.624E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.581E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    43 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  8.118E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       577.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE              14.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.118E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.111E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       584.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE             -11.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.118E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.365E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       578.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE              6.35       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.118E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.468E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       578.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE              4.39       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.118E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.666E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       584.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE             -11.4       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.118E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.522E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    43 STARTS. ATTEMPT NUMBER  3, TIME INCREMENT  2.030E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       572.       TIME AVG. FORCE        507.    
 LARGEST RESIDUAL FORCE             0.429       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.030E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.596E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION    3.742E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION
          TIME INCREMENT MAY NOW INCREASE TO    3.044E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.030E-03,  FRACTION OF STEP COMPLETED  0.369    
 STEP TIME COMPLETED       0.369    ,  TOTAL TIME COMPLETED        0.369    


  INCREMENT    44 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.044E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       575.       TIME AVG. FORCE        509.    
 LARGEST RESIDUAL FORCE            -0.510       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.044E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.174E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.567E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.044E-03,  FRACTION OF STEP COMPLETED  0.372    
 STEP TIME COMPLETED       0.372    ,  TOTAL TIME COMPLETED        0.372    


  INCREMENT    45 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.567E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       580.       TIME AVG. FORCE        510.    
 LARGEST RESIDUAL FORCE            -0.554       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.567E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.699E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    6.850E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.567E-03,  FRACTION OF STEP COMPLETED  0.377    
 STEP TIME COMPLETED       0.377    ,  TOTAL TIME COMPLETED        0.377    


  INCREMENT    46 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.850E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       586.       TIME AVG. FORCE        512.    
 LARGEST RESIDUAL FORCE              4.04       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.850E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.403E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       592.       TIME AVG. FORCE        512.    
 LARGEST RESIDUAL FORCE             -10.7       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.850E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.368E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       586.       TIME AVG. FORCE        512.    
 LARGEST RESIDUAL FORCE              5.85       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.850E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.104E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       586.       TIME AVG. FORCE        512.    
 LARGEST RESIDUAL FORCE              3.11       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.850E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.164E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       595.       TIME AVG. FORCE        512.    
 LARGEST RESIDUAL FORCE             -17.4       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.850E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.966E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    46 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.712E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       697.       TIME AVG. FORCE        514.    
 LARGEST RESIDUAL FORCE            -2.184E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.712E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.124E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.569E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.712E-03,  FRACTION OF STEP COMPLETED  0.379    
 STEP TIME COMPLETED       0.379    ,  TOTAL TIME COMPLETED        0.379    


  INCREMENT    47 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.569E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       583.       TIME AVG. FORCE        516.    
 LARGEST RESIDUAL FORCE             0.227       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.569E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.057E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.853E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.569E-03,  FRACTION OF STEP COMPLETED  0.381    
 STEP TIME COMPLETED       0.381    ,  TOTAL TIME COMPLETED        0.381    


  INCREMENT    48 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.853E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       587.       TIME AVG. FORCE        517.    
 LARGEST RESIDUAL FORCE             0.489       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.853E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.445E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.780E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.853E-03,  FRACTION OF STEP COMPLETED  0.385    
 STEP TIME COMPLETED       0.385    ,  TOTAL TIME COMPLETED        0.385    


  INCREMENT    49 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.780E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       593.       TIME AVG. FORCE        519.    
 LARGEST RESIDUAL FORCE              4.70       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.780E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.034E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       592.       TIME AVG. FORCE        519.    
 LARGEST RESIDUAL FORCE              25.1       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.780E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.308E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       663.       TIME AVG. FORCE        520.    
 LARGEST RESIDUAL FORCE             -123.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.780E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.266E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       661.       TIME AVG. FORCE        520.    
 LARGEST RESIDUAL FORCE             -116.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.780E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.306E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       675.       TIME AVG. FORCE        520.    
 LARGEST RESIDUAL FORCE             -135.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.780E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.202E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    49 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.445E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       589.       TIME AVG. FORCE        519.    
 LARGEST RESIDUAL FORCE             9.218E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.445E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.729E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.167E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.445E-03,  FRACTION OF STEP COMPLETED  0.387    
 STEP TIME COMPLETED       0.387    ,  TOTAL TIME COMPLETED        0.387    


  INCREMENT    50 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.167E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       591.       TIME AVG. FORCE        520.    
 LARGEST RESIDUAL FORCE            -6.343E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.167E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.088E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.251E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.167E-03,  FRACTION OF STEP COMPLETED  0.389    
 STEP TIME COMPLETED       0.389    ,  TOTAL TIME COMPLETED        0.389    


  INCREMENT    51 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.251E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       594.       TIME AVG. FORCE        522.    
 LARGEST RESIDUAL FORCE             6.726E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.251E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.299E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.876E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.251E-03,  FRACTION OF STEP COMPLETED  0.392    
 STEP TIME COMPLETED       0.392    ,  TOTAL TIME COMPLETED        0.392    


  INCREMENT    52 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.876E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       602.       TIME AVG. FORCE        523.    
 LARGEST RESIDUAL FORCE             -5.98       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.911E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       601.       TIME AVG. FORCE        523.    
 LARGEST RESIDUAL FORCE             -4.56       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.493E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       601.       TIME AVG. FORCE        523.    
 LARGEST RESIDUAL FORCE             -3.81       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.567E-06   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       601.       TIME AVG. FORCE        523.    
 LARGEST RESIDUAL FORCE             -3.19       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.105E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       601.       TIME AVG. FORCE        523.    
 LARGEST RESIDUAL FORCE             -2.67       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.594E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       600.       TIME AVG. FORCE        523.    
 LARGEST RESIDUAL FORCE             -2.23       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.851E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.876E-03,  FRACTION OF STEP COMPLETED  0.397    
 STEP TIME COMPLETED       0.397    ,  TOTAL TIME COMPLETED        0.397    


  INCREMENT    53 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.876E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       651.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE             -86.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.562E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       660.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE             -99.2       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.305E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       645.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -71.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.529E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       629.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -39.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.167E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       617.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -19.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.168E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       614.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -13.4       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.553E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       613.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -11.1       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.226E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       612.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -9.74       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.892E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       611.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -8.73       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.904E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                       611.       TIME AVG. FORCE        525.    
 LARGEST RESIDUAL FORCE             -7.90       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.255E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED BASED ON ALTERNATIVE
          RESIDUAL CRITERION
 

 ***WARNING: FORCE EQUILIBRIUM ACCEPTED USING THE ALTERNATE TOLERANCE.

 ITERATION SUMMARY FOR THE INCREMENT:  10 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND 10 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.876E-03,  FRACTION OF STEP COMPLETED  0.402    
 STEP TIME COMPLETED       0.402    ,  TOTAL TIME COMPLETED        0.402    


  INCREMENT    54 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.876E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              23.9       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.284E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              15.6       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.793E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              10.4       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.861E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              6.99       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.156E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              4.73       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.087E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              3.21       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.397E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       612.       TIME AVG. FORCE        526.    
 LARGEST RESIDUAL FORCE              2.19       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.418E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.876E-03,  FRACTION OF STEP COMPLETED  0.407    
 STEP TIME COMPLETED       0.407    ,  TOTAL TIME COMPLETED        0.407    


  INCREMENT    55 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.876E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       617.       TIME AVG. FORCE        528.    
 LARGEST RESIDUAL FORCE              7.17       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.251E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       617.       TIME AVG. FORCE        528.    
 LARGEST RESIDUAL FORCE              4.11       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.922E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       617.       TIME AVG. FORCE        528.    
 LARGEST RESIDUAL FORCE              2.37       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.644E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.876E-03,  FRACTION OF STEP COMPLETED  0.411    
 STEP TIME COMPLETED       0.411    ,  TOTAL TIME COMPLETED        0.411    


  INCREMENT    56 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.876E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       622.       TIME AVG. FORCE        530.    
 LARGEST RESIDUAL FORCE              1.62       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.876E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.725E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.315E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.876E-03,  FRACTION OF STEP COMPLETED  0.416    
 STEP TIME COMPLETED       0.416    ,  TOTAL TIME COMPLETED        0.416    


  INCREMENT    57 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.315E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       628.       TIME AVG. FORCE        531.    
 LARGEST RESIDUAL FORCE             0.623       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.315E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.001E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.097E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.315E-03,  FRACTION OF STEP COMPLETED  0.424    
 STEP TIME COMPLETED       0.424    ,  TOTAL TIME COMPLETED        0.424    


  INCREMENT    58 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.097E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       638.       TIME AVG. FORCE        533.    
 LARGEST RESIDUAL FORCE             -4.29       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.097E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.161E-04   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       636.       TIME AVG. FORCE        533.    
 LARGEST RESIDUAL FORCE            -0.392       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.097E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.105E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.646E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.097E-02,  FRACTION OF STEP COMPLETED  0.435    
 STEP TIME COMPLETED       0.435    ,  TOTAL TIME COMPLETED        0.435    


  INCREMENT    59 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.646E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       645.       TIME AVG. FORCE        535.    
 LARGEST RESIDUAL FORCE              6.19       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.646E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.651E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       646.       TIME AVG. FORCE        535.    
 LARGEST RESIDUAL FORCE             -2.01       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.646E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.946E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.469E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.646E-02,  FRACTION OF STEP COMPLETED  0.451    
 STEP TIME COMPLETED       0.451    ,  TOTAL TIME COMPLETED        0.451    


  INCREMENT    60 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.469E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       654.       TIME AVG. FORCE        537.    
 LARGEST RESIDUAL FORCE              66.8       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.469E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.089E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       902.       TIME AVG. FORCE        541.    
 LARGEST RESIDUAL FORCE             -354.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.469E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.801E-03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.054E+03   TIME AVG. FORCE        544.    
 LARGEST RESIDUAL FORCE             1.570E+03   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -3.555E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.679E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.036E+03   TIME AVG. FORCE        627.    
 LARGEST RESIDUAL FORCE            -6.787E+03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.848E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.287E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    60 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  6.172E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       647.       TIME AVG. FORCE        537.    
 LARGEST RESIDUAL FORCE              4.38       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.172E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.811E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       647.       TIME AVG. FORCE        537.    
 LARGEST RESIDUAL FORCE              3.24       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.172E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.221E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       647.       TIME AVG. FORCE        537.    
 LARGEST RESIDUAL FORCE              2.39       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.172E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.132E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    9.258E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.172E-03,  FRACTION OF STEP COMPLETED  0.457    
 STEP TIME COMPLETED       0.457    ,  TOTAL TIME COMPLETED        0.457    


  INCREMENT    61 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.258E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       650.       TIME AVG. FORCE        539.    
 LARGEST RESIDUAL FORCE              2.39       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         9.258E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.234E-04   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION    2.921E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION
          TIME INCREMENT MAY NOW INCREASE TO    1.389E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.258E-03,  FRACTION OF STEP COMPLETED  0.467    
 STEP TIME COMPLETED       0.467    ,  TOTAL TIME COMPLETED        0.467    


  INCREMENT    62 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.389E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       654.       TIME AVG. FORCE        541.    
 LARGEST RESIDUAL FORCE              3.23       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.389E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.574E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       654.       TIME AVG. FORCE        541.    
 LARGEST RESIDUAL FORCE             -1.32       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.389E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.346E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.083E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.389E-02,  FRACTION OF STEP COMPLETED  0.480    
 STEP TIME COMPLETED       0.480    ,  TOTAL TIME COMPLETED        0.480    


  INCREMENT    63 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.083E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       658.       TIME AVG. FORCE        543.    
 LARGEST RESIDUAL FORCE              32.5       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.083E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.273E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       755.       TIME AVG. FORCE        544.    
 LARGEST RESIDUAL FORCE             -155.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.083E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.667E-03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       806.       TIME AVG. FORCE        545.    
 LARGEST RESIDUAL FORCE              650.       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -2.506E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.122E-03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.859E+03   TIME AVG. FORCE        593.    
 LARGEST RESIDUAL FORCE            -3.898E+03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.083E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.024E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    63 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  5.207E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       655.       TIME AVG. FORCE        543.    
 LARGEST RESIDUAL FORCE              4.90       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.207E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.210E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       655.       TIME AVG. FORCE        543.    
 LARGEST RESIDUAL FORCE              3.42       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.207E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.409E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       655.       TIME AVG. FORCE        543.    
 LARGEST RESIDUAL FORCE              2.38       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.207E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.018E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.811E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.207E-03,  FRACTION OF STEP COMPLETED  0.486    
 STEP TIME COMPLETED       0.486    ,  TOTAL TIME COMPLETED        0.486    


  INCREMENT    64 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.811E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       656.       TIME AVG. FORCE        544.    
 LARGEST RESIDUAL FORCE              5.37       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.811E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.062E-04   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       656.       TIME AVG. FORCE        544.    
 LARGEST RESIDUAL FORCE              3.73       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.811E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.900E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       656.       TIME AVG. FORCE        544.    
 LARGEST RESIDUAL FORCE              2.42       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.811E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.385E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.172E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.811E-03,  FRACTION OF STEP COMPLETED  0.493    
 STEP TIME COMPLETED       0.493    ,  TOTAL TIME COMPLETED        0.493    


  INCREMENT    65 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.172E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       658.       TIME AVG. FORCE        546.    
 LARGEST RESIDUAL FORCE              2.03       AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.172E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.594E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION    4.503E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION
          TIME INCREMENT MAY NOW INCREASE TO    1.758E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.172E-02,  FRACTION OF STEP COMPLETED  0.505    
 STEP TIME COMPLETED       0.505    ,  TOTAL TIME COMPLETED        0.505    


  INCREMENT    66 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.758E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       661.       TIME AVG. FORCE        548.    
 LARGEST RESIDUAL FORCE              22.1       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.758E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.651E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       709.       TIME AVG. FORCE        549.    
 LARGEST RESIDUAL FORCE             -78.7       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.758E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.808E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       694.       TIME AVG. FORCE        548.    
 LARGEST RESIDUAL FORCE              251.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.885E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.950E-03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.512E+03   TIME AVG. FORCE        561.    
 LARGEST RESIDUAL FORCE            -1.056E+03   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.758E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.174E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    66 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  4.394E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       658.       TIME AVG. FORCE        548.    
 LARGEST RESIDUAL FORCE              4.58       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.394E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.014E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       658.       TIME AVG. FORCE        548.    
 LARGEST RESIDUAL FORCE              2.82       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.394E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.683E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       658.       TIME AVG. FORCE        548.    
 LARGEST RESIDUAL FORCE              1.74       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.394E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.024E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    6.591E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.394E-03,  FRACTION OF STEP COMPLETED  0.510    
 STEP TIME COMPLETED       0.510    ,  TOTAL TIME COMPLETED        0.510    


  INCREMENT    67 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.591E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       666.       TIME AVG. FORCE        550.    
 LARGEST RESIDUAL FORCE             -12.4       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.591E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.736E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       659.       TIME AVG. FORCE        549.    
 LARGEST RESIDUAL FORCE              26.3       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.591E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.004E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       665.       TIME AVG. FORCE        550.    
 LARGEST RESIDUAL FORCE             -13.4       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.591E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.885E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       670.       TIME AVG. FORCE        550.    
 LARGEST RESIDUAL FORCE             -21.7       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.591E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.000E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       876.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE             1.238E+03   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.591E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.624E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    67 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.648E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       659.       TIME AVG. FORCE        549.    
 LARGEST RESIDUAL FORCE              1.44       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.648E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.037E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION    9.772E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION
          TIME INCREMENT MAY NOW INCREASE TO    2.472E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.648E-03,  FRACTION OF STEP COMPLETED  0.511    
 STEP TIME COMPLETED       0.511    ,  TOTAL TIME COMPLETED        0.511    


  INCREMENT    68 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.472E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              5.70       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.242E-06   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              5.50       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.164E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              5.29       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.701E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              5.07       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.363E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              4.86       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.117E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              4.66       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.955E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              4.48       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.816E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              4.31       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.688E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              4.14       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.573E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                       659.       TIME AVG. FORCE        551.    
 LARGEST RESIDUAL FORCE              3.99       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.468E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED BASED ON ALTERNATIVE
          RESIDUAL CRITERION
 

 ***WARNING: FORCE EQUILIBRIUM ACCEPTED USING THE ALTERNATE TOLERANCE.

 ITERATION SUMMARY FOR THE INCREMENT:  10 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND 10 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.472E-03,  FRACTION OF STEP COMPLETED  0.514    
 STEP TIME COMPLETED       0.514    ,  TOTAL TIME COMPLETED        0.514    


  INCREMENT    69 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.472E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              6.41       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.551E-06   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              6.16       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.109E-06   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              5.93       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.390E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              5.72       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.128E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              5.52       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.933E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              5.37       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.780E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              5.23       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.652E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              5.08       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.543E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              4.95       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.447E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                       659.       TIME AVG. FORCE        553.    
 LARGEST RESIDUAL FORCE              4.81       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.360E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED BASED ON ALTERNATIVE
          RESIDUAL CRITERION
 

 ***WARNING: FORCE EQUILIBRIUM ACCEPTED USING THE ALTERNATE TOLERANCE.

 ITERATION SUMMARY FOR THE INCREMENT:  10 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND 10 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.472E-03,  FRACTION OF STEP COMPLETED  0.516    
 STEP TIME COMPLETED       0.516    ,  TOTAL TIME COMPLETED        0.516    


  INCREMENT    70 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.472E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              5.06       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.664E-06   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              4.84       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.342E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              4.61       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.078E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              4.38       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.953E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              4.15       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.873E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              3.93       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.809E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              3.71       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.748E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              3.50       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.685E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              3.30       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.621E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                       659.       TIME AVG. FORCE        554.    
 LARGEST RESIDUAL FORCE              3.10       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.555E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED BASED ON ALTERNATIVE
          RESIDUAL CRITERION
 

 ***WARNING: FORCE EQUILIBRIUM ACCEPTED USING THE ALTERNATE TOLERANCE.

 ITERATION SUMMARY FOR THE INCREMENT:  10 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND 10 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.472E-03,  FRACTION OF STEP COMPLETED  0.519    
 STEP TIME COMPLETED       0.519    ,  TOTAL TIME COMPLETED        0.519    


  INCREMENT    71 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.472E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       659.       TIME AVG. FORCE        556.    
 LARGEST RESIDUAL FORCE             0.495       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.889E-06   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.472E-03,  FRACTION OF STEP COMPLETED  0.521    
 STEP TIME COMPLETED       0.521    ,  TOTAL TIME COMPLETED        0.521    


  INCREMENT    72 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.472E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       661.       TIME AVG. FORCE        557.    
 LARGEST RESIDUAL FORCE             -2.60       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.472E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.375E-06   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.707E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.472E-03,  FRACTION OF STEP COMPLETED  0.524    
 STEP TIME COMPLETED       0.524    ,  TOTAL TIME COMPLETED        0.524    


  INCREMENT    73 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.707E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       664.       TIME AVG. FORCE        559.    
 LARGEST RESIDUAL FORCE             -7.99       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.707E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.382E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       663.       TIME AVG. FORCE        559.    
 LARGEST RESIDUAL FORCE             -4.83       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.707E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.910E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       661.       TIME AVG. FORCE        558.    
 LARGEST RESIDUAL FORCE             -2.93       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.707E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.798E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       661.       TIME AVG. FORCE        558.    
 LARGEST RESIDUAL FORCE             -1.70       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.707E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.099E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.561E-03

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.707E-03,  FRACTION OF STEP COMPLETED  0.527    
 STEP TIME COMPLETED       0.527    ,  TOTAL TIME COMPLETED        0.527    


  INCREMENT    74 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.561E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       661.       TIME AVG. FORCE        560.    
 LARGEST RESIDUAL FORCE              7.26       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.561E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.420E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       661.       TIME AVG. FORCE        560.    
 LARGEST RESIDUAL FORCE              1.56       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.561E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.261E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    8.341E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.561E-03,  FRACTION OF STEP COMPLETED  0.533    
 STEP TIME COMPLETED       0.533    ,  TOTAL TIME COMPLETED        0.533    


  INCREMENT    75 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.341E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       662.       TIME AVG. FORCE        561.    
 LARGEST RESIDUAL FORCE              1.09       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.341E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.569E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.251E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.341E-03,  FRACTION OF STEP COMPLETED  0.541    
 STEP TIME COMPLETED       0.541    ,  TOTAL TIME COMPLETED        0.541    


  INCREMENT    76 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.251E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       664.       TIME AVG. FORCE        563.    
 LARGEST RESIDUAL FORCE              13.7       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.251E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.700E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       675.       TIME AVG. FORCE        563.    
 LARGEST RESIDUAL FORCE             -19.1       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.251E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.743E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       665.       TIME AVG. FORCE        563.    
 LARGEST RESIDUAL FORCE              25.0       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.251E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.689E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       685.       TIME AVG. FORCE        563.    
 LARGEST RESIDUAL FORCE             -35.5       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.251E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.076E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    76 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.128E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       663.       TIME AVG. FORCE        563.    
 LARGEST RESIDUAL FORCE              1.62       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.128E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.661E-06   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.692E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.128E-03,  FRACTION OF STEP COMPLETED  0.544    
 STEP TIME COMPLETED       0.544    ,  TOTAL TIME COMPLETED        0.544    


  INCREMENT    77 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.692E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       663.       TIME AVG. FORCE        564.    
 LARGEST RESIDUAL FORCE              2.51       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.692E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.640E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.038E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.692E-03,  FRACTION OF STEP COMPLETED  0.549    
 STEP TIME COMPLETED       0.549    ,  TOTAL TIME COMPLETED        0.549    


  INCREMENT    78 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.038E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       664.       TIME AVG. FORCE        565.    
 LARGEST RESIDUAL FORCE              2.48       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         7.038E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.749E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.056E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.038E-03,  FRACTION OF STEP COMPLETED  0.556    
 STEP TIME COMPLETED       0.556    ,  TOTAL TIME COMPLETED        0.556    


  INCREMENT    79 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.056E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       666.       TIME AVG. FORCE        566.    
 LARGEST RESIDUAL FORCE              1.28       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.056E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.400E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.584E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.056E-02,  FRACTION OF STEP COMPLETED  0.567    
 STEP TIME COMPLETED       0.567    ,  TOTAL TIME COMPLETED        0.567    


  INCREMENT    80 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.584E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       669.       TIME AVG. FORCE        568.    
 LARGEST RESIDUAL FORCE              24.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.584E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.681E-04   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       723.       TIME AVG. FORCE        568.    
 LARGEST RESIDUAL FORCE             -87.1       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.584E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.269E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       709.       TIME AVG. FORCE        568.    
 LARGEST RESIDUAL FORCE              275.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.740E-02   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.033E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.646E+03   TIME AVG. FORCE        580.    
 LARGEST RESIDUAL FORCE            -1.208E+03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.584E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.270E-02   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    80 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.959E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       666.       TIME AVG. FORCE        568.    
 LARGEST RESIDUAL FORCE              3.51       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.959E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.623E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       666.       TIME AVG. FORCE        568.    
 LARGEST RESIDUAL FORCE              2.25       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.959E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.175E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.938E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-03,  FRACTION OF STEP COMPLETED  0.571    
 STEP TIME COMPLETED       0.571    ,  TOTAL TIME COMPLETED        0.571    


  INCREMENT    81 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.938E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       675.       TIME AVG. FORCE        569.    
 LARGEST RESIDUAL FORCE             -14.9       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.938E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.401E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       683.       TIME AVG. FORCE        569.    
 LARGEST RESIDUAL FORCE             -28.1       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.938E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.805E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       677.       TIME AVG. FORCE        569.    
 LARGEST RESIDUAL FORCE             -17.9       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.938E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.233E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       668.       TIME AVG. FORCE        569.    
 LARGEST RESIDUAL FORCE             -2.50       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.938E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.824E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    8.907E-03

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.938E-03,  FRACTION OF STEP COMPLETED  0.576    
 STEP TIME COMPLETED       0.576    ,  TOTAL TIME COMPLETED        0.576    


  INCREMENT    82 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.907E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       669.       TIME AVG. FORCE        570.    
 LARGEST RESIDUAL FORCE            -0.264       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         8.907E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.760E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.336E-02

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.907E-03,  FRACTION OF STEP COMPLETED  0.585    
 STEP TIME COMPLETED       0.585    ,  TOTAL TIME COMPLETED        0.585    


  INCREMENT    83 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.336E-02
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       671.       TIME AVG. FORCE        571.    
 LARGEST RESIDUAL FORCE              15.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.336E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.594E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       689.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -31.9       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.336E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.962E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       674.       TIME AVG. FORCE        571.    
 LARGEST RESIDUAL FORCE              64.8       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.336E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.869E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       763.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -147.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.336E-02   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.758E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    83 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.340E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       668.       TIME AVG. FORCE        571.    
 LARGEST RESIDUAL FORCE              17.3       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.613E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       841.       TIME AVG. FORCE        573.    
 LARGEST RESIDUAL FORCE              934.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.116E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       926.       TIME AVG. FORCE        574.    
 LARGEST RESIDUAL FORCE             1.298E+03   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.447E-03   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.750E-03   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       970.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             -634.       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -7.481E-03   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.008E-03   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       714.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -98.7       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.946E-03   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       714.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -81.7       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.120E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       708.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -66.9       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.034E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       701.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -54.4       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.180E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       685.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -30.3       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.055E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                       690.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE             -40.5       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.340E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.412E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    83 STARTS. ATTEMPT NUMBER  3, TIME INCREMENT  1.670E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       669.       TIME AVG. FORCE        571.    
 LARGEST RESIDUAL FORCE             0.680       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.670E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.252E-05   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.505E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.670E-03,  FRACTION OF STEP COMPLETED  0.587    
 STEP TIME COMPLETED       0.587    ,  TOTAL TIME COMPLETED        0.587    


  INCREMENT    84 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.505E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       669.       TIME AVG. FORCE        572.    
 LARGEST RESIDUAL FORCE              1.08       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.505E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.024E-06   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.758E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.505E-03,  FRACTION OF STEP COMPLETED  0.590    
 STEP TIME COMPLETED       0.590    ,  TOTAL TIME COMPLETED        0.590    


  INCREMENT    85 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.758E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       670.       TIME AVG. FORCE        574.    
 LARGEST RESIDUAL FORCE              3.14       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.758E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.556E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       670.       TIME AVG. FORCE        574.    
 LARGEST RESIDUAL FORCE              1.78       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.758E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.044E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.637E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.758E-03,  FRACTION OF STEP COMPLETED  0.593    
 STEP TIME COMPLETED       0.593    ,  TOTAL TIME COMPLETED        0.593    


  INCREMENT    86 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.637E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       673.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             -9.34       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.637E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.686E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       683.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             -23.8       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.637E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.731E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       671.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE              6.38       AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.637E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.342E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       680.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             -18.2       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.637E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.507E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       678.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             -15.5       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.637E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.538E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       686.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             -28.6       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.637E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.335E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    86 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.409E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       670.       TIME AVG. FORCE        575.    
 LARGEST RESIDUAL FORCE             0.576       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.409E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.187E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.114E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.409E-03,  FRACTION OF STEP COMPLETED  0.595    
 STEP TIME COMPLETED       0.595    ,  TOTAL TIME COMPLETED        0.595    


  INCREMENT    87 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.114E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       670.       TIME AVG. FORCE        576.    
 LARGEST RESIDUAL FORCE            -0.611       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.114E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.708E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.171E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.114E-03,  FRACTION OF STEP COMPLETED  0.597    
 STEP TIME COMPLETED       0.597    ,  TOTAL TIME COMPLETED        0.597    


  INCREMENT    88 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.171E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       672.       TIME AVG. FORCE        577.    
 LARGEST RESIDUAL FORCE             -2.18       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.171E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.854E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.756E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.171E-03,  FRACTION OF STEP COMPLETED  0.600    
 STEP TIME COMPLETED       0.600    ,  TOTAL TIME COMPLETED        0.600    


  INCREMENT    89 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.756E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       671.       TIME AVG. FORCE        578.    
 LARGEST RESIDUAL FORCE              9.90       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.756E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.101E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       671.       TIME AVG. FORCE        578.    
 LARGEST RESIDUAL FORCE              23.8       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.756E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.182E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       673.       TIME AVG. FORCE        578.    
 LARGEST RESIDUAL FORCE             -5.21       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.756E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.059E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       673.       TIME AVG. FORCE        578.    
 LARGEST RESIDUAL FORCE             -4.36       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.756E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.766E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       671.       TIME AVG. FORCE        578.    
 LARGEST RESIDUAL FORCE              34.8       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.756E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.364E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    89 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.189E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       671.       TIME AVG. FORCE        578.    
 LARGEST RESIDUAL FORCE            -0.746       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.189E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.564E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION   -4.648E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION
          TIME INCREMENT MAY NOW INCREASE TO    1.784E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.189E-03,  FRACTION OF STEP COMPLETED  0.601    
 STEP TIME COMPLETED       0.601    ,  TOTAL TIME COMPLETED        0.601    


  INCREMENT    90 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.784E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       671.       TIME AVG. FORCE        579.    
 LARGEST RESIDUAL FORCE             0.785       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.784E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.658E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.675E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.784E-03,  FRACTION OF STEP COMPLETED  0.603    
 STEP TIME COMPLETED       0.603    ,  TOTAL TIME COMPLETED        0.603    


  INCREMENT    91 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.675E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       671.       TIME AVG. FORCE        580.    
 LARGEST RESIDUAL FORCE              1.93       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.675E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.222E-06   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.013E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.675E-03,  FRACTION OF STEP COMPLETED  0.606    
 STEP TIME COMPLETED       0.606    ,  TOTAL TIME COMPLETED        0.606    


  INCREMENT    92 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.013E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       672.       TIME AVG. FORCE        581.    
 LARGEST RESIDUAL FORCE              4.84       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.013E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.067E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       672.       TIME AVG. FORCE        581.    
 LARGEST RESIDUAL FORCE              3.51       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.013E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.387E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       672.       TIME AVG. FORCE        581.    
 LARGEST RESIDUAL FORCE              2.54       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         4.013E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.073E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    6.019E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.013E-03,  FRACTION OF STEP COMPLETED  0.610    
 STEP TIME COMPLETED       0.610    ,  TOTAL TIME COMPLETED        0.610    


  INCREMENT    93 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.019E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       697.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE             -47.2       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.697E-04   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       673.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE              9.87       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.825E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       685.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE             -24.4       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.651E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       673.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE              7.99       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.316E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       673.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE              6.47       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.796E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       673.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE              3.66       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.597E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       702.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE             -50.1       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.019E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.338E-06   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    93 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.505E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       672.       TIME AVG. FORCE        582.    
 LARGEST RESIDUAL FORCE             0.882       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.505E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.848E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION    6.883E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION
          TIME INCREMENT MAY NOW INCREASE TO    2.257E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.505E-03,  FRACTION OF STEP COMPLETED  0.611    
 STEP TIME COMPLETED       0.611    ,  TOTAL TIME COMPLETED        0.611    


  INCREMENT    94 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.257E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       673.       TIME AVG. FORCE        583.    
 LARGEST RESIDUAL FORCE            -0.894       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.257E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.009E-05   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.386E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.257E-03,  FRACTION OF STEP COMPLETED  0.613    
 STEP TIME COMPLETED       0.613    ,  TOTAL TIME COMPLETED        0.613    


  INCREMENT    95 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.386E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       855.       TIME AVG. FORCE        586.    
 LARGEST RESIDUAL FORCE             -397.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.855E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.171E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       804.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -236.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.387E-03   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       763.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -150.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.061E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       749.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -129.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.335E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       734.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -106.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.008E-04   AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       737.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -113.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.475E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       736.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -112.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.134E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       734.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -107.       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.451E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                       733.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -106.       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.386E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.207E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    95 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.693E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       673.       TIME AVG. FORCE        584.    
 LARGEST RESIDUAL FORCE            -0.905       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.693E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.057E-05   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.539E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.693E-03,  FRACTION OF STEP COMPLETED  0.615    
 STEP TIME COMPLETED       0.615    ,  TOTAL TIME COMPLETED        0.615    


  INCREMENT    96 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.539E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       674.       TIME AVG. FORCE        585.    
 LARGEST RESIDUAL FORCE             -1.03       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.539E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.268E-06   AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.809E-03

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.539E-03,  FRACTION OF STEP COMPLETED  0.618    
 STEP TIME COMPLETED       0.618    ,  TOTAL TIME COMPLETED        0.618    


  INCREMENT    97 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.809E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       674.       TIME AVG. FORCE        586.    
 LARGEST RESIDUAL FORCE              3.00       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.809E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.147E-05   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       674.       TIME AVG. FORCE        586.    
 LARGEST RESIDUAL FORCE              2.35       AT NODE          5   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.809E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.573E-06   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.714E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.809E-03,  FRACTION OF STEP COMPLETED  0.621    
 STEP TIME COMPLETED       0.621    ,  TOTAL TIME COMPLETED        0.621    


  INCREMENT    98 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.714E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       705.       TIME AVG. FORCE        587.    
 LARGEST RESIDUAL FORCE             -54.4       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.728E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       700.       TIME AVG. FORCE        587.    
 LARGEST RESIDUAL FORCE             -46.9       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.739E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       675.       TIME AVG. FORCE        587.    
 LARGEST RESIDUAL FORCE              10.1       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.168E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       695.       TIME AVG. FORCE        587.    
 LARGEST RESIDUAL FORCE             -36.6       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.407E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       676.       TIME AVG. FORCE        587.    
 LARGEST RESIDUAL FORCE             -10.4       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.485E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       675.       TIME AVG. FORCE        587.    
 LARGEST RESIDUAL FORCE              2.52       AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.510E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION    3.658E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.714E-03,  FRACTION OF STEP COMPLETED  0.627    
 STEP TIME COMPLETED       0.627    ,  TOTAL TIME COMPLETED        0.627    


  INCREMENT    99 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.714E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       675.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE             -2.45       AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.532E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF DISP.    CORRECTION   -1.638E-05
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.714E-03,  FRACTION OF STEP COMPLETED  0.633    
 STEP TIME COMPLETED       0.633    ,  TOTAL TIME COMPLETED        0.633    


  INCREMENT   100 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.714E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       699.       TIME AVG. FORCE        589.    
 LARGEST RESIDUAL FORCE             -42.6       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.245E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       687.       TIME AVG. FORCE        589.    
 LARGEST RESIDUAL FORCE             -22.3       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.036E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       676.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE              15.0       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.158E-05   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       713.       TIME AVG. FORCE        589.    
 LARGEST RESIDUAL FORCE             -17.2       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.263E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       676.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE              10.4       AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.369E-05   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                       676.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE              9.57       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.606E-04   AT NODE          6   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                       676.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE              4.05       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.017E-04   AT NODE          4   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                       689.       TIME AVG. FORCE        589.    
 LARGEST RESIDUAL FORCE             -26.2       AT NODE          6   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         5.714E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.511E-05   AT NODE          8   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT   100 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.428E-03
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                       675.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE             0.560       AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.428E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.073E-04   AT NODE          2   DOF  2
   INSTANCE: PART-1-1                                                                        
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT
 
	SYMMETRIC DIRECT SPARSE SOLVER RUNNING ON
 	1 HOST:        1 MPI RANK  x 1 THREAD
 	NUMBER OF EQUATIONS:     24
 	NUMBER OF RHS:           1
 	NUMBER OF FLOPS:         4.900e+03

     CHECK POINT  START OF SOLVER     

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =            0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       675.       TIME AVG. FORCE        588.    
 LARGEST RESIDUAL FORCE             0.189       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.428E-03   AT NODE          5   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.550E-06   AT NODE          2   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.143E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.428E-03,  FRACTION OF STEP COMPLETED  0.634    
 STEP TIME COMPLETED       0.634    ,  TOTAL TIME COMPLETED        0.634    
 

 ***ERROR: TOO MANY INCREMENTS NEEDED TO COMPLETE THE STEP

 ***ERROR: THE ANALYSIS HAS BEEN TERMINATED DUE TO PREVIOUS ERRORS. ALL OUTPUT 
           REQUESTS HAVE BEEN WRITTEN FOR THE LAST CONVERGED INCREMENT.



     ANALYSIS SUMMARY:
     TOTAL OF        100  INCREMENTS
                      24  CUTBACKS IN AUTOMATIC INCREMENTATION
                     423  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     423  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     423  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       0  WARNING MESSAGES DURING USER INPUT PROCESSING
                       6  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       2  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       2  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =      3.2    
       SYSTEM TIME (SEC)    =     0.40    
       TOTAL CPU TIME (SEC) =      3.6    
       WALLCLOCK TIME (SEC) =            7
       MEMORY PEAK (GB)     =            0
