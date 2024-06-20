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
!C******************************************************************************
!C  VARIABLES REQUIRED BY ABAQUS (THE ARGUMENTS OF THE SUBROUTINE)
!C  FOR A DESCRIPTION OF THE LIST OF VARIABLES SEE ABAQUS MANUAL (VOL. VI)
    subroutine vumat(jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
        stepTime, totalTime, dt, cmname, coordMp, charLength, &
        props, density, strainInc, relSpinInc, &
        tempOld, stretchOld, defgradOld, fieldOld, &
        stressOld, stateOld, enerInternOld, enerInelasOld, &
        tempNew, stretchNew, defgradNew, fieldNew, &
        stressNew, stateNew, enerInternNew, enerInelasNew )
!c     c
        include 'vaba_param.inc'
        INTEGER, PARAMETER :: PREC = 8
!c     c
        dimension jblock(*), props(nprops),density(*), coordMp(*), &
            charLength(*), strainInc(*), &
            relSpinInc(*), tempOld(*), &
            stretchOld(*), &
            defgradOld(*), &
            fieldOld(*), stressOld(*), &
            stateOld(*), enerInternOld(*), &
            enerInelasOld(*), tempNew(*), &
            stretchNew(*), &
            defgradNew(*), &
            fieldNew(*), &
            stressNew(*), stateNew(*), &
            enerInternNew(*), enerInelasNew(*)
!c     c
        character*80 cmname
     
        parameter (i_umt_nblock = 1, i_umt_npt = 2, i_umt_layer = 3, &
            i_umt_kspt = 4, i_umt_noel = 5 )
     
        call vumatXtrArg(jblock(i_umt_nblock), &
            ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
            stepTime, totalTime, dt, cmname, coordMp, charLength, &
            props, density, strainInc, relSpinInc, &
            tempOld, stretchOld, defgradOld, fieldOld, &
            stressOld, stateOld, enerInternOld, enerInelasOld, &
            tempNew, stretchNew, defgradNew, fieldNew, &
            stressNew, stateNew, enerInternNew, enerInelasNew, &
            jblock(i_umt_noel), jblock(i_umt_npt), &
            jblock(i_umt_layer), jblock(i_umt_kspt))
     
        return
    end subroutine vumat
     
!c     ----------------------------------------------------------------------------------
     
    subroutine vumatXtrArg (nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
        stepTime, totalTime, timeinc, cmname, coordMp, charLength, &
        props, density, strainInc, relSpinInc, &
        tempOld, stretchOld, defgradOld, fieldOld, &
        stressOld, stateOld, enerInternOld, enerInelasOld, &
        tempNew, stretchNew, defgradNew, fieldNew, &
        stressNew, stateNew, enerInternNew, enerInelasNew, &
        nElement, nMatPoint, nLayer, nSecPoint)
!c
        include 'vaba_param.inc'
        INTEGER, PARAMETER :: PREC = 8
     
!c      all arrays dimensioned by (*) are not used in this algorithm
        dimension props(nprops), density(nblock), &
            strainInc(nblock,ndir+nshr), &
            relSpinInc(nblock,nshr), defgradOld(nblock,9), &
            stressOld(nblock,ndir+nshr), &
            stateOld(nblock,nstatev), enerInternOld(nblock), &
            enerInelasOld(nblock), &
            stretchNew(nblock,ndir+nshr), defgradNew(nblock,9), &
            stressNew(nblock,ndir+nshr)
     
        dimension enerInelasNew(nblock),stateNew(nblock,nstatev), &
            enerInternNew(nblock),eigVal(nblock,3)
     
        dimension nElement(nblock),nMatPoint(nblock),nLayer(nblock), &
            nSecPoint(nblock)
     
        character*80 cmname
!c      CHARACTER*8 CMNAME
		
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
    REAL(PREC), PARAMETER::ZTOL=1.0e-10

!C  elastic constants
	REAL(PREC) :: EMOD, ENU
!C  COMPLIANCE TENSOR
	REAL(PREC),DIMENSION(ndir + nshr,ndir + nshr)::SCOMP, DDSDDE

!C  HARDENING PARAMETERS
	REAL(PREC)::AA,BB,CC
!C  HARDENING VALUES
	REAL(PREC):: HF, HPF, HFS

!C  STRESS TENSOR AND ITS INCREMENTS
	REAL(PREC),DIMENSION(ndir + nshr)::SIGMA, DSIGMA, D2SIGMA, STRESS, DSTRAN

!C  EQUIVALENT PLASTIC STRAIN AND ITS INCREMENTS
	REAL(PREC):: EPBAR, DEPBAR, D2EPBAR

!C  YIELD FUNCTION VALUE, GRADIENT AND HESSIAN
	REAL(PREC):: YF, YFS
	REAL(PREC),DIMENSION(ndir + nshr)::GYF
	REAL(PREC),DIMENSION(ndir + nshr,ndir + nshr)::HYF

!C  CONVERGENCE TOLERANCES
!    REAL(PREC),PARAMETER::TOL1=1.0E-006, TOL2=1.0E-008
	REAL(PREC),PARAMETER::TOL1=1.0E-004

!C  TEMPORARY HOLDERS
	REAL(PREC)::TT, TTA, TTB, ZALPHA, F1, FZERO, TDEPBAR, EBULK3, &
        EG2, EG, EG3, ELAM
	REAL(PREC),DIMENSION(ndir + nshr)::YVECTOR, F2, TDSIGMA
	REAL(PREC),DIMENSION(ndir + nshr,ndir + nshr)::XIMAT
	REAL(PREC),DIMENSION(ndir + nshr, ndir + nshr)::BV, IDENTITY
	REAL(PREC),DIMENSION(ndir + nshr)::ZZ
    REAL(PREC)::VT

!C  LOOP COUNTERS
    INTEGER::K1,K2,NRK,KK,LL,MM,II,JJ,NTENS, k, i

!C  NEWTON-RAPHSON MAXIMUM NUMBER OF ITERATIONS
    INTEGER,PARAMETER:: NRMAX=50
	
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
	NKMAT = NPROPS - 7
	
	ALLOCATE (KMATERIAL(NKMAT))
	KMATERIAL = PROPS(8:7+NKMAT)
			
!C!********************************************
!C RECOVER THE EQUIVALENT PLASTIC STRAIN AT THE BEGINING OF THE INCREMENT
    NTENS = ndir + nshr

!C**********************************************************************
!C     WRITE TO A FILE
!C**********************************************************************


      
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
    
    DO K1=1, ndir
        DO K2=1, ndir
            DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
    END DO
    DO K1=ndir+1, NTENS
        DDSDDE(K1, K1)=EG
    END DO 

    TTA=1.0D0/EMOD
    TTB=-ENU/EMOD 
    TT=TWO*(ONE+ENU)/EMOD
    SCOMP = ZERO

    DO K1=1, ndir
        DO K2=1, ndir
            SCOMP(K2, K1) = TTB
        END DO
        SCOMP(K1, K1) = TTA
    END DO
    DO K1=ndir+1, NTENS
        SCOMP(K1, K1)=TT
    END DO   

    
!C COMPUTE THE TRIAL STRESS : SIGMA_{N+1} = SIGMA_{N} + C[DELTA_EPSILON]

    DO k=1, nblock, 1
        EPBAR = stateOld(k,1)
        STRESS = stressOld(k,:)
        DSTRAN = strainInc(k,:)
        DSIGMA = ZERO
        DEPBAR = ZERO

        DO K1=1, NTENS
            TT = ZERO
            DO K2=1, NTENS
                TT = TT + DDSDDE(K1, K2) * DSTRAN(K2)
            END DO
            DSIGMA(K1)= TT
            SIGMA(K1) = STRESS(K1) + TT
        END DO
    !C  write(*,*)"DS",DSTRAN
    !C  write(*,*)"S",STRESS
    !C  write(*,*)"SIG",SIGMA

    !C      DO K1=1,NTENS,1
    !C      TT=DDSDDE(K1,1)*DSTRAN(1)+DDSDDE(K1,2)*DSTRAN(2)+DDSDDE(K1,3)*DSTRAN(3)
    !C      DSIGMA(K1) = TT
    !C      SIGMA(K1)=STRESS(K1)+TT
    !C      END DO

    !C!***********************TEST ZONE*********************
        !C  write(*,*)"TEST 1", ZERO ** 2
        !C  write(*,*)"TEST 2", ZERO ** (-1)
        !C  write(*,*)"TEST 3", ZERO * ZERO ** (-1)
        !C  write(*,*)"TEST 4", ZERO * ZERO ** (-1) + ONE
        !C  write(*,*)"TEST 5", 1.0E-10 ** (-2) * ZERO


    !C CHECK YIELDING CONDITION
        CALL KHARD(HF,HPF,EPBAR,AA,BB,CC)
        CALL YFUNCTION(SIGMA,NTENS,YF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)

    !C  ELASTIC STEP :  UPDATE STRESS
        IF (YF <= HF) THEN
            !C****write(80,1)TIME, EPBAR, SIGMA(1),YF, HF, AA, BB, CC
            stressNew(k,:) = SIGMA
            stateNew(k,1) = EPBAR
    !C  DDSDDE HAS BEEN DEFINED ABOVE
    !C  THE EQUIVALENT PLASTIC STRAIN, STATEV(1), REMAINS UNCHANGED
        ELSE

        !c  write(*,*)"Je passe en plastique"
        !c  IF (ISNAN(STRESS(1))) THEN
        !c    write(*,*)"Plantage"
        !c  END IF
        !c write(*,*)"STRESS", STRESS
        !c write(*,*)"SIGMA", SIGMA

    !C***********************************************
    !C MAIN LOOP : RETURN MAPPING ALGORITHM

    !C  DEFINE COMPLIANCE (note that it outputs ENGINEERING shears)
        

        !C  SCOMP CHECKED
            
            
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

            !C  write(*,*)"XIMAT", XIMAT

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

            !C  write(*,*)"BV", BV
            
            DO II=1, NTENS, 1
                ZZ(II) = YVECTOR(II)
                DO JJ = 1, II-1, 1
                    ZZ(II) = ZZ(II) - BV(II, JJ) * ZZ(JJ)
                END DO
                ZZ(II) = ZZ(II) / BV(II, II)
            END DO

            !C  write(*,*)"ZZ", ZZ
            
            DO II=NTENS, 1, -1
                D2SIGMA(II) = ZZ(II)
                DO JJ =II+1 , NTENS, 1
                    D2SIGMA(II) = D2SIGMA(II) - BV(JJ, II) * D2SIGMA(JJ)
                END DO
                D2SIGMA(II) = D2SIGMA(II) / BV(II, II)
            END DO

            !C  write(*,*)"D2SIGMA", D2SIGMA

            

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

            !C  write(*,*)"YF", YF
            !C  write(*,*)"GYF", GYF
            !C  write(*,*)"HYF", HYF
            
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
                !c  write(80,1)F1/TOL1  
                IF(FZERO<TOL1) EXIT

                !C  write(*,*)"ZZ PRE", NRK, ZZ
                
                !C  write(*,*)"YVECTOR PRE", NRK, XIMAT
                !C  write(*,*)"XIMAT PRE", NRK, XIMAT

        !C      ASSEMBLE XIMAT MATRIX AND Y-VECTOR
                DO K1=1,NTENS,1
                    YVECTOR(K1)=-(F1*GYF(K1)+HPF*F2(K1))
                    DO K2=K1,NTENS,1
                        TT=HPF*(SCOMP(K1,K2)+DEPBAR*HYF(K1,K2))+GYF(K1)*GYF(K2)
                        XIMAT(K1,K2)=TT
                        XIMAT(K2,K1)=TT
                    END DO
                END DO

                !C  write(*,*)"XIMAT", NRK, XIMAT

                DO K1=1, NTENS, 1
                    DO K2=1, NTENS, 1
                        BV(K1,K2)=XIMAT(K1,K2)
                    END DO
                END DO

                !C  write(*,*)"BV PRE", NRK, BV
                
        !C      SOLVE FOR STRESS NR-INCREMENT USING CHOLESKY ALGORITHM
                DO JJ=1, NTENS, 1
                    !C  write(*,*)"JJ, NRK, BV(JJ,JJ) PRE BOUCLE", JJ, NRK, BV(JJ, JJ)
                    DO KK=1, JJ-1, 1
                        !C  write(*,*)"JJ,KK", JJ, KK
                        BV(JJ,JJ)= BV(JJ,JJ) - BV(JJ,KK) * BV(JJ,KK)
                    END DO
                    !C  write(*,*)"JJ, NRK, BV(JJ,JJ) POST BOUCLE", JJ, NRK, BV(JJ, JJ)
                    BV(JJ,JJ) = DSQRT(BV(JJ,JJ))
                    !C  write(*,*)"JJ, NRK, BV(JJ,JJ) POST BOUCLE ET SQRT", JJ, NRK, BV(JJ, JJ)
                    DO II=(JJ+1), NTENS, 1
                        DO KK=1, JJ-1, 1
                            BV(II,JJ) = BV(II,JJ) - BV(II,KK) * BV(JJ,KK)
                        END DO
                        BV(II,JJ)=BV(II,JJ)/BV(JJ,JJ)
                    END DO
                END DO
                
                !C  write(*,*)"BV POST", NRK, BV

                DO II=1, NTENS, 1
                    ZZ(II) = YVECTOR(II)
                    DO JJ = 1, II-1, 1
                        ZZ(II) = ZZ(II) - BV(II, JJ) * ZZ(JJ)
                    END DO
                    ZZ(II) = ZZ(II) / BV(II, II)
                END DO
                
                !C  write(*,*)"ZZ POST", NRK, ZZ
                !C  write(*,*)"BV POST", NRK, BV
                !C  write(*,*)"YVECTOR POST", NRK, XIMAT
                !C  write(*,*)"XIMAT POST", NRK, XIMAT


                DO II=NTENS, 1, -1
                    D2SIGMA(II) = ZZ(II)
                    DO JJ =II+1 , NTENS, 1
                        D2SIGMA(II) = D2SIGMA(II) - BV(JJ, II) * D2SIGMA(JJ)
                    END DO
                    D2SIGMA(II) = D2SIGMA(II) / BV(II, II)
                END DO

        !C      CALCULATE EQUIVALENT PLASTIC STRAIN NR-INCREMENT 
                D2EPBAR=F1
                !C  write(*,*),"D2PBAR PRE", NRK, D2EPBAR
                !C  write(*,*)"GYF", NRK, GYF
                !C  write(*,*)"D2SIGMA POST", NRK, D2SIGMA
                !C  write(*,*)"HPF", NRK, HPF


                DO K1=1,NTENS,1
                    D2EPBAR=D2EPBAR+GYF(K1)*D2SIGMA(K1)
                END DO
                D2EPBAR=D2EPBAR/HPF

        !C      DO LINE SEARCH
                TDEPBAR=DEPBAR+D2EPBAR
                TDSIGMA=DSIGMA+D2SIGMA

                !C  write(*,*),"D2PBAR POST", NRK, D2EPBAR

                CALL LSEARCH(NTENS,STRESS,TDSIGMA,DSTRAN,EPBAR,TDEPBAR,FZERO, &
                        SCOMP,KMATERIAL,NKMAT,AA,BB,CC,ZALPHA,DEGREE,NCOEFF,NMON)

        !C      UPDATE

                DEPBAR=DEPBAR+ZALPHA*D2EPBAR
                DSIGMA=DSIGMA+ZALPHA*D2SIGMA
                
                !C  write(*,*),"DEPBAR", NRK, DEPBAR

            END DO !!! END OF NEWTON-RAPHSON ITERATIONS
                
        !C  UPDATE STATE VARIABLE

            !C  write(*,*)"EPBAR, DEPBAR", EPBAR, DEPBAR
            stateNew(k,1)=EPBAR+DEPBAR

        !C  UPDATE STRESS
            stressNew(k,:) = STRESS+DSIGMA
        !C  write(*,*)"DSIGMA", DSIGMA
    !C************************************** COMPUTE TANGENT MODULUS: DDSDDE
        END IF
    
    END DO
    do i=1, nblock
		stress_power = ONE/TWO * ( &
            (stressOld(i,1) + stressNew(i,1)) * strainInc(i,1) + &
            (stressOld(i,2) + stressNew(i,2)) * strainInc(i,2) + &
            (stressOld(i,3) + stressNew(i,3)) * strainInc(i,3) + &
            TWO * (stressOld(i,4) + stressNew(i,4)) * strainInc(i,4) + &
            TWO * (stressOld(i,5) + stressNew(i,5)) * strainInc(i,5) + &
            TWO * (stressOld(i,6) + stressNew(i,6)) * strainInc(i,6) )
        enerInternNew(i) = enerInternOld(i) + stress_power/density(i)
            
        smean = ONE/THREE * &
            (stressNew(i,1) + stressNew(i,2) + stressNew(i,3))
        equiv_stress = SQRT(THREE/TWO * &
                ((stressNew(i,1) - smean) ** 2 + &
                (stressNew(i,2) - smean) ** 2 + &
                (stressNew(i,3) - smean) ** 2 + &
                TWO * stressNew(i,4) ** 2 + &
                TWO * stressNew(i,5) ** 2 + &
                TWO * stressNew(i,6) ** 2))
        plastic_work_inc = equiv_stress * (stateNew(1, i)-stateOld(1, i))
        enerInelasNew(i) = enerInelasOld(i) + plastic_work_inc / density(i)
	END DO
        !C  write(*,*)"DDSDDE", DDSDDE
    DEALLOCATE(KMATERIAL)	
    RETURN
    END SUBROUTINE vumatXtrArg


!C**************************HARDENING***************************
!C*****NOTE:
!C*****THIS UMAT IDENTIFIES THE HARDENING SUB BY THE NAME 'KHARD'
!C*****(DEACTIVATE THE OTHER BY RENAMING)
     	
!C****: Swift (Power hardening law)
    SUBROUTINE   KHARD(HF,HPF,EPBAR,AAZ,BBZ,CCZ)
!C      COMPUTES THE HARDENING AND ITS DERIVATIVE
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: PREC = 8
    REAL(PREC) :: HF, HPF, EPBAR, AAZ, BBZ, CCZ
    HF  = AAZ*((BBZ+EPBAR)**CCZ)
    HPF =  (CCZ/(BBZ+EPBAR))*HF
    RETURN
    END SUBROUTINE  KHARD

!C****: Voce (Exponential hardening law)
    SUBROUTINE  swKHARD(HF,HPF,EPBAR,AAZ,BBZ,CCZ)
!C      COMPUTES THE HARDENING AND ITS DERIVATIVE
    
    IMPLICIT NONE
    INTEGER, PARAMETER :: PREC = 8
    REAL(PREC) :: HF, HPF, EPBAR, AAZ, BBZ, CCZ, ONE
    ONE = 1.0D0
    HF= BBZ * (ONE - EXP(-CCZ*EPBAR))
	HPF = -CCZ * BBZ * EXP(-CCZ * EPBAR)
    HF  = AAZ-HF
    RETURN
    END SUBROUTINE  swKHARD

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
    REAL(PREC),DIMENSION(NTENS - 1)::DEVIA
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL
	REAL(PREC)::YF, MAX
	REAL(PREC)::ZTOL=1.0E-10
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
    DEVIA(3)=SIGMA(4)
    DEVIA(4)=SIGMA(5)
    DEVIA(5)=SIGMA(6)



	YF = 0.0D0
	
    N0 = 1
	DO MM=0,DEGREE,1
        DO LL=0, DEGREE - MM, 1
            DO KK=0, DEGREE - MM - LL, 1
                DO JJ=0, DEGREE - MM - LL - KK, 1
                    II = DEGREE - MM - LL - KK - JJ
                    IF (((KK==LL) .AND. (LL==MM)) .OR. & 
                    ((MOD(KK,2)==0) .AND. (MOD(LL,2)==0) &
                        .AND. (MOD(MM,2)==0))) THEN

                        YF = YF + KMATERIAL(N0) &
                        * DEVIA(1) ** DBLE(II) &
                        * DEVIA(2) ** DBLE(JJ) &
                        * DEVIA(3) ** DBLE(KK) &
                        * DEVIA(4) ** DBLE(LL) &
                        * DEVIA(5) ** DBLE(MM)
                        N0 = N0 + 1
                    END IF
                END DO
            END DO
        END DO
    END DO

    YF = YF**(1.0D0/DBLE(DEGREE)) 

	RETURN
	END SUBROUTINE YFUNCTION

    SUBROUTINE GYFUNCTION(SIGMA,NTENS,YF,GYF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

    INTEGER::NTENS,NKMAT,DEGREE,NCOEFF,NMON
	REAL(PREC),DIMENSION(NTENS)::SIGMA, GYF
	REAL(PREC),DIMENSION(NTENS - 1)::GYFP
    REAL(PREC),DIMENSION(NTENS - 1)::DEVIA
	REAL(PREC),DIMENSION(NTENS - 1, NTENS)::GDEV
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL
	REAL(PREC)::YF,ZYF,MAX
    REAL(PREC)::ZTOL=1.0E-10
    
    REAL(PREC), PARAMETER::ZERO=0.0D0
    REAL(PREC), PARAMETER::ONE=1.0D0
    REAL(PREC), PARAMETER::TWO=2.0D0
    REAL(PREC), PARAMETER::THREE=3.0D0
    REAL(PREC), PARAMETER::SIX=6.0D0
    
	INTEGER::II,JJ,KK,LL,MM,K1,K2,N0

	
    DEVIA(1)=SIGMA(1) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(2)=SIGMA(2) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(3)=SIGMA(4)
    DEVIA(4)=SIGMA(5)
    DEVIA(5)=SIGMA(6)


	
    YF = 0.0D0
	GYF = 0.0D0
	GYFP = 0.0D0
	GDEV = 0.0D0
	
	GDEV(1,1) = TWO / THREE
	GDEV(1,2) = - ONE / THREE
	GDEV(1,3) = - ONE / THREE
	GDEV(2,1) = - ONE / THREE
	GDEV(2,2) = TWO / THREE
	GDEV(2,3) = - ONE / THREE
	GDEV(3,4) = ONE
	GDEV(4,5) = ONE
	GDEV(5,6) = ONE
	
    !c  write(*,*)"DEVIA", DEVIA
    N0 = 1
    DO MM=0,DEGREE,1
        DO LL=0, DEGREE - MM, 1
            DO KK=0, DEGREE - LL - MM, 1
                DO JJ=0, DEGREE - LL - MM - KK, 1
                    II=DEGREE - LL - MM - KK - JJ
                    IF (((KK==LL) .AND. (LL==MM)) .OR. & 
                    ((MOD(KK,2)==0) .AND. (MOD(LL,2)==0) &
                    .AND. (MOD(MM,2)==0))) THEN
                        YF = YF + KMATERIAL(N0) &
                            * DEVIA(1) ** DBLE(II) &
                            * DEVIA(2) ** DBLE(JJ) &
                            * DEVIA(3) ** DBLE(KK) &
                            * DEVIA(4) ** DBLE(LL) &
                            * DEVIA(5) ** DBLE(MM)
							
							
							
						IF (II > 0) THEN
							GYFP(1) = GYFP(1) + KMATERIAL(N0) &
								* DBLE(II) * DEVIA(1) ** (DBLE(II-1)) &
								* DEVIA(2) ** DBLE(JJ) &
								* DEVIA(3) ** DBLE(KK) &
								* DEVIA(4) ** DBLE(LL) &
								* DEVIA(5) ** DBLE(MM)
						END IF
						IF (JJ > 0) THEN
							GYFP(2) = GYFP(2) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DBLE(JJ) * DEVIA(2) ** (DBLE(JJ-1)) &
								* DEVIA(3) ** DBLE(KK) &
								* DEVIA(4) ** DBLE(LL) &
								* DEVIA(5) ** DBLE(MM)
						END IF
						IF (KK > 0) THEN
							GYFP(3) = GYFP(3) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DEVIA(2) ** DBLE(JJ) &
								* DBLE(KK) * DEVIA(3) ** (DBLE(KK-1)) &
								* DEVIA(4) ** DBLE(LL) &
								* DEVIA(5) ** DBLE(MM)
						END IF
						IF (LL > 0) THEN
							GYFP(4) = GYFP(4) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DEVIA(2) ** DBLE(JJ) &
								* DEVIA(3) ** DBLE(KK) &
								* DBLE(LL) * DEVIA(4) ** (DBLE(LL-1)) &
								* DEVIA(5) ** DBLE(MM)
						END IF
						IF (MM > 0) THEN
							GYFP(5) = GYFP(5) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DEVIA(2) ** DBLE(JJ) &
								* DEVIA(3) ** DBLE(KK) &
								* DEVIA(4) ** DBLE(LL) &
								* DBLE(MM) * DEVIA(5) ** (DBLE(MM-1))
						END IF
                        !C write(*,*)"DEVIA", DEVIA
                        !C write(*,*)"DBLE(II), JJ, KK, LL, MM", DBLE(II), JJ, KK, LL, MM
                        !C write(*,*)"DEVIA(1) ** DBLE(II)", DEVIA(1) ** DBLE(II)
                        !C write(*,*)"DEVIA(2) ** JJ", DEVIA(2) ** JJ
                        !C write(*,*)"DEVIA(4) ** KK", DEVIA(4) ** KK
                        !C write(*,*)"DEVIA(5) ** LL", DEVIA(5) ** LL
                        !C write(*,*)"DEVIA(6) ** (MM-1)", DEVIA(6) ** (MM-1)
                        !C write(*,*)"GYF(6)", GYF(6)
                        N0 = N0 + 1
                    END IF
                END DO
            END DO
        END DO
    END DO
	
	DO K1=1, NTENS
		DO K2= 1, NTENS - 1
			GYF(K1) = GYF(K1) + GYFP(K2) * GDEV(K2, K1)
		END DO
	END DO
	
    ZYF = YF ** (ONE/DBLE(DEGREE))

    DO II = 1, NTENS
        GYF(II) = GYF(II) * (ZYF/(DBLE(DEGREE) * YF))
    END DO
    
    YF = ZYF

	RETURN
    
	END SUBROUTINE GYFUNCTION
!CC***************************************************************************

    SUBROUTINE HYFUNCTION(SIGMA,NTENS,YF,GYF,HYF,KMATERIAL,NKMAT,DEGREE,NCOEFF,NMON)
    IMPLICIT NONE
	INTEGER,PARAMETER::PREC=8

	INTEGER::NTENS,NKMAT,DEGREE,NCOEFF,NMON
	REAL(PREC),DIMENSION(NTENS)::SIGMA,GYF
	REAL(PREC),DIMENSION(NTENS - 1)::GYFP
	REAL(PREC),DIMENSION(NTENS - 1)::DEVIA
	REAL(PREC),DIMENSION(NTENS - 1, NTENS)::GDEV
	REAL(PREC),DIMENSION(NTENS,NTENS)::HYF
	REAL(PREC),DIMENSION(NTENS -1,NTENS -1)::HYFP
	REAL(PREC),DIMENSION(NKMAT)::KMATERIAL
	REAL(PREC)::YF, MAX
    REAL(PREC)::ZTOL=1.0E-10
	REAL(PREC)::ZYF,YVAL,Y2VAL
	INTEGER::II,JJ,KK,MM,LL,K1,K2,K3,K4,N0
    REAL(PREC), PARAMETER::ZERO=0.0D0
    REAL(PREC), PARAMETER::ONE=1.0D0
    REAL(PREC), PARAMETER::TWO=2.0D0
    REAL(PREC), PARAMETER::THREE=3.0D0
    REAL(PREC), PARAMETER::SIX=6.0D0


	DEVIA(1) = SIGMA(1) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(2) = SIGMA(2) - ONE/THREE * (SIGMA(1) + SIGMA(2)&
                + SIGMA(3))
    DEVIA(3) = SIGMA(4)
    DEVIA(4) = SIGMA(5)
    DEVIA(5) = SIGMA(6)
    
    YF = 0.0D0
	GYF = 0.0D0
	GYFP = 0.0D0
    HYF = 0.0D0
	HYFP = 0.0D0
	
	GDEV(1,1) = TWO / THREE
	GDEV(1,2) = - ONE / THREE
	GDEV(1,3) = - ONE / THREE
	GDEV(2,1) = - ONE / THREE
	GDEV(2,2) = TWO / THREE
	GDEV(2,3) = - ONE / THREE
	GDEV(3,4) = ONE
	GDEV(4,5) = ONE
	GDEV(5,6) = ONE
    
    N0 = 1
    DO MM=0,DEGREE,1
        DO LL=0, DEGREE - MM, 1
            DO KK=0, DEGREE - LL - MM, 1
                DO JJ=0, DEGREE - LL - MM - KK, 1
                    II=DEGREE - LL - MM - KK - JJ
                    IF (((KK==LL) .AND. (LL==MM)) .OR. & 
                    ((MOD(KK,2)==0) .AND. (MOD(LL,2)==0) &
                    .AND. (MOD(MM,2)==0))) THEN
						
                        YF = YF + KMATERIAL(N0) &
                            * DEVIA(1) ** DBLE(II) &
                            * DEVIA(2) ** DBLE(JJ) &
                            * DEVIA(3) ** DBLE(KK) &
                            * DEVIA(4) ** DBLE(LL) &
                            * DEVIA(5) ** DBLE(MM)
							
							
							
							
							
							
						IF (II > 0) THEN
							GYFP(1) = GYFP(1) + KMATERIAL(N0) &
								* DBLE(II) * DEVIA(1) ** (DBLE(II-1)) &
								* DEVIA(2) ** DBLE(JJ) &
								* DEVIA(3) ** DBLE(KK) &
								* DEVIA(4) ** DBLE(LL) &
								* DEVIA(5) ** DBLE(MM)
							IF (II > 1) THEN
								HYFP(1,1) = HYFP(1,1) + KMATERIAL(N0) &
									* DBLE(II) * (DBLE(II-1)) * DEVIA(1) ** (DBLE(II-2)) &
									* DEVIA(2) ** DBLE(JJ) &
									* DEVIA(3) ** DBLE(KK) &
									* DEVIA(4) ** DBLE(LL) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (JJ > 0) THEN
								HYFP(1,2) = HYFP(1,2) + KMATERIAL(N0) &
									* DBLE(II) * DEVIA(1) ** (DBLE(II-1)) &
									* DBLE(JJ) * DEVIA(2) ** (DBLE(JJ-1)) &
									* DEVIA(3) ** DBLE(KK) &
									* DEVIA(4) ** DBLE(LL) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (KK > 0) THEN
								HYFP(1,3) = HYFP(1,3) + KMATERIAL(N0) &
									* DBLE(II) * DEVIA(1) ** (DBLE(II-1)) &
									* DEVIA(2) ** DBLE(JJ) &
									* DBLE(KK) * DEVIA(3) ** (DBLE(KK - 1)) &
									* DEVIA(4) ** DBLE(LL) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (LL > 0) THEN
								HYFP(1,4) = HYFP(1,4) + KMATERIAL(N0) &
									* DBLE(II) * DEVIA(1) ** (DBLE(II-1)) &
									* DEVIA(2) ** DBLE(JJ) &
									* DEVIA(3) ** DBLE(KK) &
									* DBLE(LL) * DEVIA(4) ** (DBLE(LL - 1)) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (MM > 0) THEN 
								HYFP(1,5) = HYFP(1,5) + KMATERIAL(N0) &
									* DBLE(II) * DEVIA(1) ** (DBLE(II-1)) &
									* DEVIA(2) ** DBLE(JJ) &
									* DEVIA(3) ** DBLE(KK) &
									* DEVIA(4) ** DBLE(LL) &
									* DBLE(MM) * DEVIA(5) ** (DBLE(MM - 1))
							END IF
						END IF
						
						
						
						
						IF (JJ > 0) THEN
							GYFP(2) = GYFP(2) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DBLE(JJ) * DEVIA(2) ** (DBLE(JJ-1)) &
								* DEVIA(3) ** DBLE(KK) &
								* DEVIA(4) ** DBLE(LL) &
								* DEVIA(5) ** DBLE(MM)
							IF (JJ > 1) THEN
								HYFP(2,2) = HYFP(2,2) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DBLE(JJ) * (DBLE(JJ-1)) * DEVIA(2) ** (DBLE(JJ - 2)) &
									* DEVIA(3) ** DBLE(KK) &
									* DEVIA(4) ** DBLE(LL) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (KK > 0) THEN
								HYFP(2,3) = HYFP(2,3) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DBLE(JJ) * DEVIA(2) ** (DBLE(JJ - 1)) &
									* DBLE(KK) * DEVIA(3) ** (DBLE(KK - 1)) &
									* DEVIA(4) ** DBLE(LL) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (LL > 0) THEN
								HYFP(2,4) = HYFP(2,4) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DBLE(JJ) * DEVIA(2) ** (DBLE(JJ - 1)) &
									* DEVIA(3) ** DBLE(KK) &
									* DBLE(LL) * DEVIA(4) ** (DBLE(LL - 1)) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (MM > 0) THEN 
								HYFP(2,5) = HYFP(2,5) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DBLE(JJ) * DEVIA(2) ** (DBLE(JJ - 1)) &
									* DEVIA(3) ** DBLE(KK) &
									* DEVIA(4) ** DBLE(LL) &
									* DBLE(MM) * DEVIA(5) ** (DBLE(MM - 1))
							END IF
						END IF
						
						
						
						
						
						
						IF (KK > 0) THEN
							GYFP(3) = GYFP(3) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DEVIA(2) ** DBLE(JJ) &
								* DBLE(KK) * DEVIA(3) ** (DBLE(KK-1)) &
								* DEVIA(4) ** DBLE(LL) &
								* DEVIA(5) ** DBLE(MM)
							IF (KK > 1) THEN
								HYFP(3,3) = HYFP(3,3) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DEVIA(2) ** DBLE(JJ) &
									* DBLE(KK) * (DBLE(KK - 1)) * DEVIA(3) ** (DBLE(KK - 2)) &
									* DEVIA(4) ** DBLE(LL) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (LL > 0) THEN
								HYFP(3,4) = HYFP(3,4) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DEVIA(2) ** DBLE(JJ) &
									* DBLE(KK) * DEVIA(3) ** (DBLE(KK - 1)) &
									* DBLE(LL) * DEVIA(4) ** (DBLE(LL - 1)) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (MM > 0) THEN 
								HYFP(3,5) = HYFP(3,5) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DEVIA(2) ** DBLE(JJ)&
									* DBLE(KK) * DEVIA(3) ** (DBLE(KK - 1)) &
									* DEVIA(4) ** DBLE(LL) &
									* DBLE(MM) * DEVIA(5) ** (DBLE(MM - 1))
							END IF								
						END IF
						
						
						
						
						
						IF (LL > 0) THEN
							GYFP(4) = GYFP(4) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DEVIA(2) ** DBLE(JJ) &
								* DEVIA(3) ** DBLE(KK) &
								* DBLE(LL) * DEVIA(4) ** (DBLE(LL-1)) &
								* DEVIA(5) ** DBLE(MM)
							IF (LL > 1) THEN
								HYFP(4,4) = HYFP(4,4) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DEVIA(2) ** DBLE(JJ) &
									* DEVIA(3) ** DBLE(KK) &
									* DBLE(LL) * (DBLE(LL - 1)) * DEVIA(4) ** (DBLE(LL - 2)) &
									* DEVIA(5) ** DBLE(MM)
							END IF
							IF (MM > 0) THEN 
								HYFP(4,5) = HYFP(4,5) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DEVIA(2) ** DBLE(JJ) &
									* DEVIA(3) ** DBLE(KK) &
									* DBLE(LL) * DEVIA(4) ** (DBLE(LL - 1)) &
									* DBLE(MM) * DEVIA(5) ** (DBLE(MM - 1))
							END IF	
						END IF
						
						
						IF (MM > 0) THEN
							GYFP(5) = GYFP(5) + KMATERIAL(N0) &
								* DEVIA(1) ** DBLE(II) &
								* DEVIA(2) ** DBLE(JJ) &
								* DEVIA(3) ** DBLE(KK) &
								* DEVIA(4) ** DBLE(LL) &
								* DBLE(MM) * DEVIA(5) ** (DBLE(MM-1))
							IF (MM > 1) THEN 
								HYFP(5,5) = HYFP(5,5) + KMATERIAL(N0) &
									* DEVIA(1) ** DBLE(II) &
									* DEVIA(2) ** DBLE(JJ) &
									* DEVIA(3) ** DBLE(KK) &
									* DEVIA(4) ** DBLE(LL) &
									* DBLE(MM) * (DBLE(MM - 1)) * DEVIA(5) ** (DBLE(MM - 2))
							END IF	
						END IF
                        N0 = N0 + 1
                    END IF
                END DO
            END DO
        END DO
	END DO
	
	DO K1=1, NTENS
		DO K2= 1, NTENS - 1
			GYF(K1) = GYF(K1) + GYFP(K2) * GDEV(K2, K1)
		END DO
	END DO
	
	DO K1=1, NTENS
		DO K2 = 1, NTENS
			DO K3 = 1, NTENS - 1
				DO K4 = 1, NTENS - 1
					HYF(K1, K2) = HYF(K1, K2) + GDEV(K4, K1) * HYFP(K4, K3) * GDEV(K3, K2)
				END DO
			END DO
		END DO
	END DO
	
    ZYF = YF**(1.0D0/DBLE(DEGREE))
	YVAL = ZYF/(DBLE(DEGREE)*YF)
    Y2VAL = DBLE(DEGREE-1)/ZYF

    !C  write(*,*)"YF", YF
    !C  write(*,*)"GYF PRE", GYF

    !C TO CHECK AGAIN WITH FORMULA OF COMPOSED FUNCTION
    DO II = 1, NTENS
        GYF(II) = GYF(II) * (ZYF/(DBLE(DEGREE) * YF))
    END DO

    !C  write(*,*)"GYF POST", GYF
    !C  write(*,*)"Y2VAL", Y2VAL
    !C  write(*,*)"HYF PRE", HYF
    DO II = 1, NTENS
        DO JJ = II, NTENS
			HYF(II, JJ) =HYF(II, JJ) * YVAL - GYF(II) * GYF(JJ) * Y2VAL
            HYF(JJ, II) = HYF(II, JJ)
        END DO
    END DO

    !C  write(*,*)"HYF POST", HYF
    YF = ZYF 

	RETURN
    
	END SUBROUTINE HYFUNCTION

!CC**************************** END OF PolyN IMPLEMENTATION