!DIR$ FREEFORM
!CC*******************************************************************************
!C  UMAT: developed by S.C. Soare. 
!C  IT USES ISOTROPIC ELASTICITY COUPLED WITH AN ANISOTROPIC YIELD FUNCTION. 
!C  THE SUBROUTINE IS DESIGNED FOR PLANE (2D) STRESS STATES. 
!C  IT IS BASED ON THE FULLY IMPLICIT RETURN MAPPING ALGORITHM (with quadratic line search) 
!C  ONE STATE VARIABLE : HARDENING PARAMETER (THE EQUIVALENT PLASTIC STRAIN).
!C  IT IS ASSUMED THAT THE USER HAS DEFINED (USING THE *ORIENTATION OPTION IN ABAQUS)
!C  A LOCAL COORDINATE SYSTEM THAT ROTATES WITH THE MATERIAL (SO THAT DEROT IS ALWAYS
!C  EQUAL WITH THE IDENTITY MATRIX).
!C*****************************************************************************
!C The vector PROPS(NPROPS) contains the material properties as defined in the 
!C  *MATERIAL=USER option in ABAQUS in the following order
!C  PROPS(1) = EMOD  (Elasticity: Young modulus)
!C  PROPS(2) = MU  (Elasticity: Poisson ratio)
!C  Hardening laws defined: 
!C  Swift (power-)law: sigma^bar = a*(b + ep^bar)**c
!C  Voce (exp-)law: sigma^bar = a - b*exp(-c*ep^bar)  (default)
!C  Read further naming/renaming convention in the HARDENING section of this code 
!C  (more specific hardening laws can be implemented in the same section) 
!C  PROPS(3) = a 
!C  PROPS(4) = b
!C  PROPS(5) = c
!C  PROPS(6),...,PROPS(NPROPS): PARAMETERS OF YIELD FUNCTION
!C  Yield func implemented: PolyN
!C  Note: The first two parameters are the degree and number of coefficients
!C!************************************************************************************** 

!CCCC-----NOTE: the UMAT interface may vary with the ABAQUS version 	
!C	SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
!C      RPL,DDSDDT,DRPLDE,DRPLDT, &
!C      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
!C      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROTT,PNEWDT, &
!C      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

!CCC---NOTE: the INCLUDE directive enforces implicit casting in conflict with 'IMPLICIT NONE'
!CC          use 'IMPLICIT NONE' in the testing/implementation phase and then comment it out 
    IMPLICIT NONE      	
!	INCLUDE 'ABA_PARAM.INC'
!C    INTEGER, PARAMETER :: PREC =  SELECTED_REAL_KIND(15,307)
    INTEGER, PARAMETER :: PREC = 8
!C******************************************************************************
!C  VARIABLES REQUIRED BY ABAQUS (THE ARGUMENTS OF THE SUBROUTINE)
!C  FOR A DESCRIPTION OF THE LIST OF VARIABLES SEE ABAQUS MANUAL (VOL. VI)

!C    !!!CHARACTER(80)::  CMNAME
	CHARACTER*8 CMNAME
    REAL(PREC)::SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,PNEWDT,CELENT
    INTEGER::NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KINC
    REAL(PREC),DIMENSION(NTENS):: STRESS,DDSDDT,DRPLDE,STRAN,DSTRAN
    REAL(PREC),DIMENSION (NTENS, NTENS) :: DDSDDE 
    REAL(PREC),DIMENSION(NSTATV) :: STATEV
    REAL(PREC),DIMENSION(NPROPS) :: PROPS
    REAL(PREC),DIMENSION(3,3) :: DFGRD0, DFGRD1, DROT
    REAL(PREC),DIMENSION(3) :: COORDS
    REAL(PREC),DIMENSION(2) :: TIME
    REAL(PREC),DIMENSION(1) :: PREDEF, DPRED
    INTEGER,DIMENSION(4)::JSTEP

		
!C    CHARACTER*80 CMNAME
!C	INTEGER::NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KINC
!C	REAL(PREC)::SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,PNEWDT,CELENT
!C    REAL(PREC),DIMENSION::STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS), &
!C	DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS), &
!C	TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS), &
!C	COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!C	INTEGER,DIMENSION::JSTEP(4)	
	
!C*******************************************************************************
!C  INTERNAL VARIABLES OF THE SUBROUTINE
      
!C  VECTOR O MATERIAL PARAMETERS
    INTEGER::NKMAT
	REAL(PREC),ALLOCATABLE, DIMENSION(:)::KMATERIAL
    
	REAL(PREC), PARAMETER::ZERO=0.0D0
	REAL(PREC), PARAMETER::ONE=1.0D0
	REAL(PREC), PARAMETER::TWO=2.0D0
	REAL(PREC), PARAMETER::THREE=3.0D0
	REAL(PREC), PARAMETER::SIX=6.0D0

!C  elastic constants
	REAL(PREC) :: EMOD, ENU
!C  COMPLIANCE TENSOR
	REAL(PREC),DIMENSION(NTENS,NTENS)::SCOMP

!C  HARDENING PARAMETERS
	REAL(PREC)::AA,BB,CC
!C  HARDENING VALUES
	REAL(PREC):: HF, HPF

!C  STRESS TENSOR AND ITS INCREMENTS
	REAL(PREC),DIMENSION(NTENS)::SIGMA, DSIGMA, D2SIGMA

!C  EQUIVALENT PLASTIC STRAIN AND ITS INCREMENTS
	REAL(PREC):: EPBAR, DEPBAR, D2EPBAR

!C  YIELD FUNCTION VALUE, GRADIENT AND HESSIAN
	REAL(PREC):: YF
	REAL(PREC),DIMENSION(NTENS)::GYF
	REAL(PREC),DIMENSION(NTENS,NTENS)::HYF

!C  CONVERGENCE TOLERANCES
!    REAL(PREC),PARAMETER::TOL1=1.0E-006, TOL2=1.0E-008
	REAL(PREC),PARAMETER::TOL1=1.0E-007

!C  TEMPORARY HOLDERS
	REAL(PREC)::TT, TTA, TTB, ZALPHA, F1, FZERO, TDEPBAR, EBULK3, &
        EG2, EG, EG3, ELAM
	REAL(PREC),DIMENSION(NTENS)::YVECTOR, F2, TDSIGMA
	REAL(PREC),DIMENSION(NTENS,NTENS)::XIMAT
	REAL(PREC),DIMENSION(NTENS, NTENS)::BV, IDENTITY
	REAL(PREC),DIMENSION(NTENS)::ZZ

!C  LOOP COUNTERS
    INTEGER::K1,K2,NRK,KK,LL,MM,II,JJ

!C  NEWTON-RAPHSON MAXIMUM NUMBER OF ITERATIONS
    INTEGER,PARAMETER:: NRMAX=100  
	
!C PolyN interface variables 
    INTEGER::DEGREE,NCOEFF,NMON	

!C*****************************************************
	EMOD = PROPS(1)
    ENU = PROPS(2)
    AA = PROPS(3)
    BB = PROPS(4)
    CC = PROPS(5)
	DEGREE = INT(PROPS(6))
	NCOEFF = INT(PROPS(7))
	NKMAT = (NPROPS - 7)/2
    NMON=INT(DEGREE/2)
	
	ALLOCATE (KMATERIAL(NKMAT))
	KMATERIAL = PROPS(8:7+NKMAT)
			
!C!********************************************
!C RECOVER THE EQUIVALENT PLASTIC STRAIN AT THE BEGINING OF THE INCREMENT
    EPBAR = STATEV(1)
      
!C!********************************************
!C INITIALIZE THE STIFFNESS TENSOR (IT WILL BE STORED IN DDSDDE)
    DDSDDE = ZERO

!C ELASTIC PROPERTIES (3D STRESS)
    EBULK3 = EMOD/(ONE - TWO*ENU)
    EG2 = EMOD/(ONE+ENU)
    EG = EG2/TWO
    EG3 = THREE*EG
    ELAM = (EBULK3 - EG2)/THREE

!C COMPUTE THE STIFFNESS TENSOR IN 3D
    
    DO K1=1, NDI
        DO K2=1, NDI
            DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
    END DO
    DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
    END DO    
    
!C COMPUTE THE TRIAL STRESS : SIGMA_{N+1} = SIGMA_{N} + C[DELTA_EPSILON]

    DO K1=1, NTENS
        TT = ZERO
        DO K2=1, NTENS
            TT = TT + DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
        DSIGMA(K1)= TT
        SIGMA(K1) = STRESS(K1) + TT
    END DO
    write(*,*)"DS",DSTRAN
    write(*,*)"S",STRESS
    write(*,*)"SIG",SIGMA

!C      DO K1=1,NTENS,1
!C      TT=DDSDDE(K1,1)*DSTRAN(1)+DDSDDE(K1,2)*DSTRAN(2)+DDSDDE(K1,3)*DSTRAN(3)
!C      DSIGMA(K1) = TT
!C      SIGMA(K1)=STRESS(K1)+TT
!C      END DO

!C CHECK YIELDING CONDITION
    CALL KHARD(HF,HPF,EPBAR,AA,BB,CC)
    CALL YFUNCTION(SIGMA,NTENS,YF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
    
    write(*,*)"HF",HF
    write(*,*)"YF",YF

!C  ELASTIC STEP :  UPDATE STRESS
	IF (YF <= HF) THEN
	    STRESS = SIGMA
!C  DDSDDE HAS BEEN DEFINED ABOVE
!C  THE EQUIVALENT PLASTIC STRAIN, STATEV(1), REMAINS UNCHANGED
        DEALLOCATE(KMATERIAL)		
        RETURN
	END IF

!C***********************************************
!C MAIN LOOP : RETURN MAPPING ALGORITHM

!C  DEFINE COMPLIANCE (note that it outputs ENGINEERING shears)
    
    TTA=1.0D0/EMOD
	TTB=-ENU/EMOD 
    TT=TWO*(ONE+ENU)/EMOD

    DO K1=1, NDI
        DO K2=1, NDI
            SCOMP(K2, K1) = TTB
        END DO
        SCOMP(K1, K1) = TTA
    END DO
    DO K1=NDI+1, NTENS
        SCOMP(K1, K1)=TT
    END DO
    
	
!C  SCOMP(1,1)=TTA
!C  SCOMP(2,1)=TTB
!C  SCOMP(3,1)=ZERO
!C  SCOMP(1,2)=TTB
!C  SCOMP(2,2)=TTA
!C  SCOMP(3,2)=ZERO
!C  SCOMP(1,3)=ZERO
!C  SCOMP(2,3)=ZERO
!C  SCOMP(3,3)=2.0D0*(1.0D0+ENU)/EMOD

!C    FIRST Newton-Raphson step (no Hessian required)
!C**************************************************************      
!C  DEPBAR=ZERO

	CALL GYFUNCTION(SIGMA,NTENS,YF,GYF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
	F1=YF-HF

!C  ASSEMBLE XIMAT MATRIX AND Y-VECTOR
    DO K1=1,NTENS,1
	    YVECTOR(K1)=-F1*GYF(K1)
	    DO K2=K1,NTENS,1
	        TT=HPF*SCOMP(K1,K2)+GYF(K1)*GYF(K2)
	        XIMAT(K1,K2)=TT
	        XIMAT(K2,K1)=TT
	    END DO
	END DO

    DO K1=1, NTENS, 1
        DO K2=1, NTENS, 1
            BV(K1,K2)=XIMAT(K1,K2)
        END DO
    END DO

!C  SOLVE FOR STRESS NR-INCREMENT USING CHOLESKY ALGORITHM
    DO JJ=1, NTENS, 1
        DO KK=1, JJ-1, 1
            BV(JJ,JJ)= BV(JJ,JJ) - BV(JJ,KK) * BV(JJ,KK)
        END DO
        BV(JJ,JJ) = DSQRT(BV(JJ,JJ))
        DO II=(JJ+1), NTENS, 1
            DO KK=1, JJ-1, 1
                BV(II,JJ)=BV(II,JJ) - BV(II,KK) * BV(JJ,KK)
            END DO
            BV(II,JJ)=BV(II,JJ)/BV(JJ,JJ)
        END DO
    END DO
    
	DO II=1, NTENS, 1
        ZZ(II) = YVECTOR(II)
        DO JJ = 1, II-1, 1
            ZZ(II) = ZZ(II) - BV(II, JJ) * ZZ(JJ)
        END DO
        ZZ(II) = ZZ(II) / BV(II, II)
    END DO
    
    DO II=NTENS, 1, -1
        D2SIGMA(II) = ZZ(II)
        DO JJ =II+1 , NTENS, 1
            D2SIGMA(II) = D2SIGMA(II) - BV(JJ, II) * D2SIGMA(JJ)
        END DO
        D2SIGMA(II) = D2SIGMA(II) / BV(II, II)
    END DO

!C  CALCULATE EQUIVALENT PLASTIC STRAIN NR-INCREMENT 
    D2EPBAR=F1
	DO K1=1,NTENS,1
	    D2EPBAR=D2EPBAR+GYF(K1)*D2SIGMA(K1)
	END DO
	D2EPBAR=D2EPBAR/HPF
    
!C  DO LINE SEARCH (along the full NR-step)
    TDEPBAR=D2EPBAR
	TDSIGMA=DSIGMA+D2SIGMA	
    FZERO=F1
    CALL LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO, &
	             SCOMP,KMATERIAL,NKMAT,AA,BB,CC,ZALPHA,DEGREE,NCOEFF,NMON)

!C  UPDATE
    DEPBAR=ZALPHA*D2EPBAR
	DSIGMA=DSIGMA+ZALPHA*D2SIGMA
	
!C    THE REST OF N-R ITERATIONS
!C******************************************************	     
    DO NRK=1,NRMAX,1

!C      CALCULATE NEW VALUES ASSOCIATED WITH NEW STATE
        CALL KHARD(HF,HPF,EPBAR+DEPBAR,AA,BB,CC)
	    SIGMA=STRESS+DSIGMA
	    CALL HYFUNCTION(SIGMA,NTENS,YF,GYF,HYF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
		
	    F1=YF-HF
	    FZERO=F1*F1
	    DO K1=1,NTENS,1
	        TT=DEPBAR*GYF(K1)-DSTRAN(K1)
	        DO K2=1,NTENS,1
	            TT=TT+SCOMP(K1,K2)*DSIGMA(K2)
	        END DO
	        F2(K1)=TT
	        FZERO=FZERO+TT*TT
	    END DO
        FZERO=DSQRT(FZERO)
!C      CHECK TOLERANCES
!C        IF ((DABS(F1)<TOL1).AND.(DSQRT(TTB)<TOL2)) EXIT
        IF(FZERO<TOL1) EXIT

!C      ASSEMBLE XIMAT MATRIX AND Y-VECTOR
        DO K1=1,NTENS,1
	        YVECTOR(K1)=-(F1*GYF(K1)+HPF*F2(K1))
	        DO K2=K1,NTENS,1
	            TT=HPF*(SCOMP(K1,K2)+DEPBAR*HYF(K1,K2))+GYF(K1)*GYF(K2)
	            XIMAT(K1,K2)=TT
                XIMAT(K2,K1)=TT
            END DO
	    END DO

!C      SOLVE FOR STRESS NR-INCREMENT USING CHOLESKY ALGORITHM
        DO JJ=1, NTENS, 1
            DO KK=1, JJ-1, 1
                BV(JJ,JJ)= BV(JJ,JJ) - BV(JJ,KK) * BV(JJ,KK)
            END DO
            BV(JJ,JJ) = DSQRT(BV(JJ,JJ))
            DO II=(JJ+1), NTENS, 1
                DO KK=1, JJ-1, 1
                    BV(II,JJ)=BV(II,JJ) - BV(II,KK) * BV(JJ,KK)
                END DO
                BV(II,JJ)=BV(II,JJ)/BV(JJ,JJ)
            END DO
        END DO
        
        DO II=1, NTENS, 1
            ZZ(II) = YVECTOR(II)
            DO JJ = 1, II-1, 1
                ZZ(II) = ZZ(II) - BV(II, JJ) * ZZ(JJ)
            END DO
            ZZ(II) = ZZ(II) / BV(II, II)
        END DO
        
        DO II=NTENS, 1, -1
            D2SIGMA(II) = ZZ(II)
            DO JJ =II+1 , NTENS, 1
                D2SIGMA(II) = D2SIGMA(II) - BV(JJ, II) * D2SIGMA(JJ)
            END DO
            D2SIGMA(II) = D2SIGMA(II) / BV(II, II)
        END DO

!C      CALCULATE EQUIVALENT PLASTIC STRAIN NR-INCREMENT 
        D2EPBAR=F1
	    DO K1=1,NTENS,1
	        D2EPBAR=D2EPBAR+GYF(K1)*D2SIGMA(K1)
	    END DO
	    D2EPBAR=D2EPBAR/HPF

!C      DO LINE SEARCH
        TDEPBAR=DEPBAR+D2EPBAR
	    TDSIGMA=DSIGMA+D2SIGMA
        CALL LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO, &
	             SCOMP,KMATERIAL,NKMAT,AA,BB,CC,ZALPHA,DEGREE,NCOEFF,NMON)

!C      UPDATE
        DEPBAR=DEPBAR+ZALPHA*D2EPBAR
	    DSIGMA=DSIGMA+ZALPHA*D2SIGMA

	END DO !!! END OF NEWTON-RAPHSON ITERATIONS
        
!C  UPDATE STATE VARIABLE
    STATEV(1)=EPBAR+DEPBAR

!C  UPDATE STRESS
    STRESS = STRESS+DSIGMA

!C************************************** COMPUTE TANGENT MODULUS: DDSDDE

!C  COMPUTE XIMAT MATRIX 
    DO K1=1,NTENS,1
	    DO K2=K1,NTENS,1
	        TT=SCOMP(K1,K2)+DEPBAR*HYF(K1,K2)
	        XIMAT(K1,K2)=TT
	        XIMAT(K2,K1)=TT
	    END DO
	END DO

!C  INVERT XIMAT AND STORE XIMAT^(-1) INTO SCOMP (NO LONGER NEEDED)
    DO JJ=1, NTENS, 1
        DO KK=1, JJ-1, 1
            BV(JJ,JJ)= BV(JJ,JJ) - BV(JJ,KK) * BV(JJ,KK)
        END DO
        BV(JJ,JJ) = DSQRT(BV(JJ,JJ))
        DO II=(JJ+1), NTENS, 1
            DO KK=1, JJ-1, 1
                BV(II,JJ)=BV(II,JJ) - BV(II,KK) * BV(JJ,KK)
            END DO
            BV(II,JJ)=BV(II,JJ)/BV(JJ,JJ)
        END DO
    END DO
    
   IDENTITY=ZERO
    
    DO II=1, NTENS, 1
        IDENTITY(II,II)=ONE
    END DO
	
    DO LL=1, NTENS, 1
        DO II=1, NTENS, 1
            ZZ(II) = IDENTITY(II, LL)
            DO JJ = 1, II-1, 1
                ZZ(II) = ZZ(II) - BV(II, JJ) * ZZ(JJ)
            END DO
            ZZ(II) = ZZ(II) / BV(II, II)
        END DO
        
        DO II=NTENS, 1, -1
            SCOMP(II, LL) = ZZ(II)
            DO JJ =II+1 , NTENS, 1
                SCOMP(II, LL) = SCOMP(II, LL) - BV(JJ, II) * SCOMP(JJ, LL)
            END DO
            SCOMP(II, LL) = SCOMP(II, LL) / BV(II, II)
        END DO
    END DO
    
    SCOMP(1,2)=SCOMP(2,1)
    SCOMP(1,3)=SCOMP(3,1)
    SCOMP(1,4)=SCOMP(4,1)
    SCOMP(1,5)=SCOMP(5,1)
    SCOMP(1,6)=SCOMP(6,1)
    SCOMP(2,3)=SCOMP(3,2)
    SCOMP(2,4)=SCOMP(4,2)
    SCOMP(2,5)=SCOMP(5,2)
    SCOMP(2,6)=SCOMP(6,2)
    SCOMP(3,4)=SCOMP(4,3)
    SCOMP(3,5)=SCOMP(5,3)
    SCOMP(3,6)=SCOMP(6,3)
    SCOMP(4,5)=SCOMP(5,4)
    SCOMP(4,6)=SCOMP(6,4)
    SCOMP(5,6)=SCOMP(6,5)
	
	
!C  CALCULATE  SCOMP[GYF] AND STORE IT INTO DSIGMA
!C  DSIGMA=(/ZERO,ZERO,ZERO/)
	DSIGMA=ZERO
    DO K1=1,NTENS,1
	    DO K2=1,NTENS,1
	        DSIGMA(K2)=DSIGMA(K2)+SCOMP(K2,K1)*GYF(K1)
	    END DO
	END DO

!C  CALCULATE 1/K
    TT=HPF
	DO K1=1,NTENS,1
	    TT=TT+GYF(K1)*DSIGMA(K1)
	END DO

!C  UPDATE DDSDDE
    DO K1=1,NTENS,1
	    DO K2=K1,NTENS,1
	        TTB=SCOMP(K1,K2)-DSIGMA(K1)*DSIGMA(K2)/TT
	        DDSDDE(K1,K2)=TTB
	        DDSDDE(K2,K1)=TTB
	    END DO
	END DO
    
	DO K1=1,NTENS,1
	    DDSDDE(K1,K1)=SCOMP(K1,K1)-DSIGMA(K1)*DSIGMA(K1)/TT
	END DO

    DEALLOCATE(KMATERIAL)	
    RETURN
    END SUBROUTINE  UMAT


!C**************************HARDENING***************************
!C*****NOTE:
!C*****THIS UMAT IDENTIFIES THE HARDENING SUB BY THE NAME 'KHARD'
!C*****(DEACTIVATE THE OTHER BY RENAMING)
     	
!C****: Swift (Power hardening law)
    SUBROUTINE   swKHARD(HF,HPF,EPBAR,AAZ,BBZ,CCZ)
!C      COMPUTES THE HARDENING AND ITS DERIVATIVE
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: PREC = 8
    REAL(PREC) :: HF, HPF, EPBAR, AAZ, BBZ, CCZ
    HF  = AAZ*((BBZ+EPBAR)**CCZ)
    HPF =  (CCZ/(BBZ+EPBAR))*HF
    RETURN
    END SUBROUTINE  swKHARD

!C****: Voce (Exponential hardening law)
    SUBROUTINE  KHARD(HF,HPF,EPBAR,AAZ,BBZ,CCZ)
!C      COMPUTES THE HARDENING AND ITS DERIVATIVE
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: PREC = 8
    REAL(PREC) :: HF, HPF, EPBAR, AAZ, BBZ, CCZ
    HF=BBZ/EXP(CCZ*EPBAR)
	HPF = CCZ*HF
    HF  = AAZ-HF
    RETURN
    END SUBROUTINE  KHARD

!C**********************************************************
	SUBROUTINE LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO, &
	                   SCOMP,KMATERIAL,NKMAT,AAZ,BBZ,CCZ,ZALPHA,DEGREE,NCOEFF,NMON) 

	IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

    INTEGER::NTENS, NKMAT,DEGREE,NCOEFF,NMON
	REAL(PREC),DIMENSION(NTENS)::STRESS,TDSIGMA,DSTRAN
	REAL(PREC),DIMENSION(NTENS,NTENS)::SCOMP
	REAL(PREC)::EPBAR,TDEPBAR,FZERO,AAZ,BBZ,CCZ,ZALPHA
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL

!C     INTERNAL VARIABLES
    REAL(PREC),DIMENSION(NTENS)::TSIGMA,GYF  
	REAL(PREC)::HF,HPF,TEPBAR,YF,TT,FONE
	INTEGER::KK,JJ

    TSIGMA=STRESS+TDSIGMA
	TEPBAR=EPBAR+TDEPBAR
	
	CALL KHARD(HF,HPF,TEPBAR,AAZ,BBZ,CCZ)
	CALL GYFUNCTION(TSIGMA,NTENS,YF,GYF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
	
    FONE=(YF-HF)
	FONE=FONE*FONE
	DO KK=1,NTENS,1
	    TT=TDEPBAR*GYF(KK)-DSTRAN(KK)
	    DO JJ=1,NTENS,1
	        TT=TT+SCOMP(KK,JJ)*TDSIGMA(JJ)
	    END DO
	FONE=FONE+TT*TT
	END DO
	FONE=DSQRT(FONE)
	ZALPHA=1.0D0
	IF(FONE<=0.5D0*FZERO) RETURN	
	
!!	ZALPHA=0.75D0*(FZERO/(2.0D0*FONE))
    ZALPHA=0.375D0*(FZERO/FONE)	
    RETURN
    END SUBROUTINE LSEARCH	
!C*********************************************************************************

!CCCCCCC********** YIELD FUNCTION CALCULATIONS ***********************************
!CCCCCCC NOTE:   YFUNCTION RETURNS JUST YIELD FUNCTION VALUE
!CCCCCCC        GYFUNCTION RETURNS YIELD FUNCTION VALUE AND GRADIENT
!CCCCCCC        HYFUNCTION RETURNS YIELD FUNCTION VALUE, GRADIENT AND HESSIAN
!C
!CCCCCC************** PolyN YIELD FUNCTION: 
!CCCCCC************** 1) ORTHOTROPIC SYMMETRY 
!CCCCCC************** 2) HOMOGENEITY DEGREE = N (any degree)
    SUBROUTINE YFUNCTION(SIGMA,NTENS,YF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

    INTEGER::NTENS,NKMAT,DEGREE,NCOEFF,NMON
	REAL(PREC),DIMENSION(NTENS)::SIGMA
    REAL(PREC),DIMENSION(NTENS)::DEVIA
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL
	REAL(PREC)::YF
	REAL(PREC),PARAMETER::ZTOL=1.0E-007
	REAL(PREC)::BB
	INTEGER::II,JJ,KK,MM,LL,N0
    
    REAL(PREC), PARAMETER::ZERO=0.0D0
    REAL(PREC), PARAMETER::ONE=1.0D0
    REAL(PREC), PARAMETER::TWO=2.0D0
    REAL(PREC), PARAMETER::THREE=3.0D0
    REAL(PREC), PARAMETER::SIX=6.0D0
	
    
    DEVIA(1)=SIGMA(1) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(2)=SIGMA(2) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(3)=SIGMA(3) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(4)=SIGMA(4)
    DEVIA(5)=SIGMA(5)
    DEVIA(6)=SIGMA(6)
    
    
    
	YF = 0.0D0
	
    N0 = 0
	DO MM=0,DEGREE,1
        DO LL=0, DEGREE - MM, 1
            DO KK=0, DEGREE - MM - LL, 1
                DO JJ=0, DEGREE - MM - LL - KK, 1
                    II = DEGREE - MM - LL - KK - JJ
                    IF (((KK==LL) .AND. (LL==MM)) .OR. & 
                    ((MOD(KK,2)==0) .AND. (MOD(LL,2)==0) &
                        .AND. (MOD(MM,2)==0))) THEN
                        YF = YF + KMATERIAL(N0) &
                        * DEVIA(1) ** II &
                        * DEVIA(2) ** JJ &
                        * DEVIA(4) ** KK &
                        * DEVIA(5) ** LL &
                        * DEVIA(6) ** MM
                        N0 = N0 + 1
                    END IF
                END DO
            END DO
        END DO
    END DO
    
	RETURN
	END SUBROUTINE YFUNCTION

    SUBROUTINE GYFUNCTION(SIGMA,NTENS,YF,GYF, KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

    INTEGER::NTENS,NKMAT,DEGREE,NCOEFF,NMON
	REAL(PREC),DIMENSION(NTENS)::SIGMA, GYF
    REAL(PREC),DIMENSION(NTENS)::DEVIA
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL
	REAL(PREC)::YF,ZZTT,ZZRHO,ZZGAMM,MSX,MSY
	REAL(PREC),PARAMETER::ZTOL=1.0E-007
	REAL(PREC)::BB
    
    REAL(PREC), PARAMETER::ZERO=0.0D0
    REAL(PREC), PARAMETER::ONE=1.0D0
    REAL(PREC), PARAMETER::TWO=2.0D0
    REAL(PREC), PARAMETER::THREE=3.0D0
    REAL(PREC), PARAMETER::SIX=6.0D0
    
	INTEGER::II,JJ,KK,MM,LL,N0
	
    DEVIA(1)=SIGMA(1) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(2)=SIGMA(2) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(3)=SIGMA(3) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(4)=SIGMA(4)
    DEVIA(5)=SIGMA(5)
    DEVIA(6)=SIGMA(6)
    
    YF = 0.0D0
	GYF = 0.0D0
    
    DO MM=0,DEGREE,1
        DO LL=0, DEGREE - MM, 1
            DO KK=0, DEGREE - LL - MM, 1
                DO JJ=0, DEGREE - LL - MM - KK, 1
                    II=DEGREE - LL - MM - KK - JJ
                    IF (((KK==LL) .AND. (LL==MM)) .OR. & 
                    ((MOD(KK,2)==0) .AND. (MOD(LL,2)==0) &
                    .AND. (MOD(MM,2)==0))) THEN
                        YF = YF + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        GYF(1) = GYF(1) + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        GYF(2) = GYF(2) + KMATERIAL(N0) &
                            * (-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        GYF(3) = GYF(3) + KMATERIAL(N0) &
                            * (-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        GYF(4) = GYF(4) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        GYF(5) = GYF(5) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM
                        GYF(6) = GYF(6) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1)
                        N0 = N0 + 1
                    END IF
                END DO
            END DO
        END DO
    END DO
	RETURN
    
	END SUBROUTINE GYFUNCTION
!CC***************************************************************************

    SUBROUTINE HYFUNCTION(SIGMA,NTENS,YF,GYF,HYF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

	INTEGER::NTENS,NKMAT,DEGREE,NCOEFF,NMON
	REAL(PREC),DIMENSION(NTENS)::SIGMA,GYF,DEVIA
	REAL(PREC),DIMENSION(NTENS,NTENS)::HYF
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL
	REAL(PREC)::YF,ZZTT,ZZRHO,ZZGAMM,MSX,MSY
    REAL(PREC),PARAMETER::ZTOL=1.0E-007
	REAL(PREC)::ZYF,YVAL,Y2VAL,ATT,TMP,BB,D1BB,D2BB,D11BB,D22BB,D12BB
	INTEGER::II,JJ,KK,MM,LL,N0
    REAL(PREC),DIMENSION(NCOEFF)::VD1,VD2,VD11,VD22,VD12
!    REAL(PREC),DIMENSION(NMON)::VD3,VD13,VD23,VD33
	REAL(PREC),DIMENSION(NMON)::VD3(0:NMON-1),VD13(0:NMON-1),VD23(0:NMON-1),VD33(0:NMON-1)
	
    REAL(PREC), PARAMETER::ZERO=0.0D0
    REAL(PREC), PARAMETER::ONE=1.0D0
    REAL(PREC), PARAMETER::TWO=2.0D0
    REAL(PREC), PARAMETER::THREE=3.0D0
    REAL(PREC), PARAMETER::SIX=6.0D0
	
	DEVIA(1)= SIGMA(1) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(2)=SIGMA(2) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(3)=SIGMA(3) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(4)=SIGMA(4)
    DEVIA(5)=SIGMA(5)
    DEVIA(6)=SIGMA(6)
    YF = 0.0D0
	GYF = 0.0D0
    HYF=0.0D0
    
    DO MM=0,DEGREE,1
        DO LL=0, DEGREE - MM, 1
            DO KK=0, DEGREE - LL - MM, 1
                DO JJ=0, DEGREE - LL - MM - KK, 1
                    II=DEGREE - LL - MM - KK - JJ
                    IF (((KK==LL) .AND. (LL==MM)) .OR. & 
                    ((MOD(KK,2)==0) .AND. (MOD(LL,2)==0) &
                    .AND. (MOD(MM,2)==0))) THEN
                        YF = YF + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(3) ** KK &
                            * DEVIA(4) ** LL &
                            * DEVIA(5) ** MM
                        GYF(1) = GYF(1) + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        GYF(2) = GYF(2) + KMATERIAL(N0) &
                            * (-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        GYF(3) = GYF(3) + KMATERIAL(N0) &
                            * (-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        GYF(4) = GYF(4) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        GYF(5) = GYF(5) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM
                        GYF(6) = GYF(6) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1)
                        HYF(1,1) = HYF(1,1) + KMATERIAL(N0) &
                            *(TWO/THREE)* (TWO/THREE) * II * (II-1) * DEVIA(1) ** (II-2) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * (-ONE/THREE) * JJ * (JJ-1) * DEVIA(2) ** (JJ-2) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + TWO * KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(2,2) = HYF(2,2) + KMATERIAL(N0) &
                            *(-ONE/THREE)* (-ONE/THREE) * II * (II-1) * DEVIA(1) ** (II-2) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * (TWO/THREE) * JJ * (JJ-1) * DEVIA(2) ** (JJ-2) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + TWO * KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(3,3) = HYF(3,3) + KMATERIAL(N0) &
                            *(-ONE/THREE)* (-ONE/THREE) * II * (II-1) * DEVIA(1) ** (II-2) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * (-ONE/THREE) * JJ * (JJ-1) * DEVIA(2) ** (JJ-2) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + TWO * KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(1,2) = HYF(1,2) + KMATERIAL(N0) &
                            *(TWO/THREE)* (-ONE/THREE) * II * (II-1) * DEVIA(1) ** (II-2) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * (TWO/THREE) * JJ * (JJ-1) * DEVIA(2) ** (JJ-2) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(1,3) = HYF(1,3) + KMATERIAL(N0) &
                            *(TWO/THREE)* (-ONE/THREE) * II * (II-1) * DEVIA(1) ** (II-2) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * (TWO/THREE) * JJ * (JJ-1) * DEVIA(2) ** (JJ-2) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(2,3) = HYF(2,3) + KMATERIAL(N0) &
                            *(-ONE/THREE)* (-ONE/THREE) * II * (II-1) * DEVIA(1) ** (II-2) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * (-ONE/THREE) * JJ * (JJ-1) * DEVIA(2) ** (JJ-2) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(1,4) = HYF(1,4) + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM 
                        HYF(1,5) = HYF(1,5) + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM 
                        HYF(1,6) = HYF(1,6) + KMATERIAL(N0) &
                            *(TWO/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1) &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1)
                        HYF(2,4) = HYF(2,4) + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** MM &
                            * DEVIA(6) ** LL &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(2,5) = HYF(2,5) + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM 
                        HYF(2,6) = HYF(2,6) + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1) &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (TWO/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1)
                        HYF(3,4) = HYF(3,4) + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** MM &
                            * DEVIA(6) ** LL &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(3,5) = HYF(3,5) + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM 
                        HYF(3,6) = HYF(2,6) + KMATERIAL(N0) &
                            *(-ONE/THREE) * II * DEVIA(1) ** (II-1) &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1) &
                            + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * (-ONE/THREE) * JJ * DEVIA(2) ** (JJ-1) &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1)
                        HYF(4,4) = HYF(4,4) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * KK * (KK-1) * DEVIA(4) ** (KK-2) &
                            * DEVIA(5) ** LL &
                            * DEVIA(6) ** MM
                        HYF(5,5) = HYF(5,5) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * (LL-1) * DEVIA(5) ** (LL-2) &
                            * DEVIA(6) ** MM
                        HYF(6,6) = HYF(6,6) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * DEVIA(5) ** LL &
                            * MM * (MM-1) * DEVIA(6) ** (MM-2)
                        HYF(4,5) = HYF(4,5) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * LL * DEVIA(5) ** (LL-1) &
                            * DEVIA(6) ** MM
                        HYF(4,6) = HYF(4,6) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * KK * DEVIA(4) ** (KK-1) &
                            * DEVIA(5) ** LL &
                            * MM * DEVIA(6) ** (MM-1)
                        HYF(5,6) = HYF(5,6) + KMATERIAL(N0) &
                            * DEVIA(1) ** II &
                            * DEVIA(2) ** JJ &
                            * DEVIA(4) ** KK &
                            * LL * DEVIA(5) ** (LL-1) &
                            * MM * DEVIA(6) ** (MM-1)
                        HYF(2,1) = HYF(1,2)
                        HYF(3,1) = HYF(1,3)
                        HYF(3,2) = HYF(2,3)
                        HYF(4,1) = HYF(1,4)
                        HYF(5,1) = HYF(1,5)
                        HYF(6,1) = HYF(1,6)
                        HYF(4,2) = HYF(2,4)
                        HYF(5,2) = HYF(2,5)
                        HYF(6,2) = HYF(2,6)
                        HYF(4,3) = HYF(3,4)
                        HYF(5,3) = HYF(3,5)
                        HYF(6,3) = HYF(3,6)
                        HYF(5,4) = HYF(4,5)
                        HYF(6,4) = HYF(4,6)
                        HYF(6,5) = HYF(5,6)
                        N0 = N0 + 1
                    END IF
                END DO
            END DO
        END DO
    END DO
	RETURN
    
	END SUBROUTINE HYFUNCTION

!CC**************************** END OF PolyN IMPLEMENTATION
