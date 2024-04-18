

c      SUBROUTINE VUMAT(
cC READ only -
c     1 NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
c     2 STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
c     3 PROPS, DENSITY, STRAININC, RELSPININC,
c     4 TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
c     5 STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
c     6 TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
cC WRITE only -
c     7 STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
cC
c      INCLUDE 'VABA_PARAM.INC'
cC
c      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK),
c     1 CHARLENGTH(NBLOCK), STRAININC(NBLOCK, NDIR+NSHR),
c     2 RELSPININC(NBLOCK, NSHR), TEMPOLD(NBLOCK),
c     3 STRETCHOLD(NBLOCK, NDIR+NSHR),DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
c     4 FIELDOLD(NBLOCK, NFIELDV), STRESSOLD(NBLOCK, NDIR+NSHR),
c     5 STATEOLD(NBLOCK, NSTATEV), ENERINTERNOLD(NBLOCK),
c     6 ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
c     7 STRETCHNEW(NBLOCK, NDIR+NSHR),DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
c     8 FIELDNEW(NBLOCK, NFIELDV), STRESSNEW(NBLOCK,NDIR+NSHR),
c     9 STATENEW(NBLOCK, NSTATEV), ENERINTERNNEW(NBLOCK),
c     1 ENERINELASNEW(NBLOCK)
C

      subroutine vumat (
c    c Read only -
     1 jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
c    c Write only -
     7 stressNew, stateNew, enerInternNew, enerInelasNew )
c     c
      include 'vaba_param.inc'
c     c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     1 charLength(*), strainInc(*),
     2 relSpinInc(*), tempOld(*),
     3 stretchOld(*),
     4 defgradOld(*),
     5 fieldOld(*), stressOld(*),
     6 stateOld(*), enerInternOld(*),
     7 enerInelasOld(*), tempNew(*),
     8 stretchNew(*),
     9 defgradNew(*), 
     1 fieldNew(*),
     2 stressNew(*), stateNew(*),
     3 enerInternNew(*), enerInelasNew(*)
c     c
      character*80 cmname
     
      parameter (
     1 i_umt_nblock = 1,
     2 i_umt_npt = 2,
     3 i_umt_layer = 3,
     4 i_umt_kspt = 4,
     5 i_umt_noel = 5 )
     
       call vumatXtrArg ( jblock(i_umt_nblock),
     1 ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, dt, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     5 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
     7 stressNew, stateNew, enerInternNew, enerInelasNew,
     8 jblock(i_umt_noel), jblock(i_umt_npt),
     9 jblock(i_umt_layer), jblock(i_umt_kspt))
     
      return
      end
     
c     ----------------------------------------------------------------------------------
     
      subroutine vumatXtrArg (
c    c read only -
     1 nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2 stepTime, totalTime, timeinc, cmname, coordMp, charLength,
     3 props, density, strainInc, relSpinInc,
     4 tempOld, stretchOld, defgradOld, fieldOld,
     3 stressOld, stateOld, enerInternOld, enerInelasOld,
     6 tempNew, stretchNew, defgradNew, fieldNew,
c    c write only -
     5 stressNew, stateNew, enerInternNew, enerInelasNew,
c    c read only extra arguments -
     6 nElement, nMatPoint, nLayer, nSecPoint)
c
       include 'vaba_param.inc'
     
c      all arrays dimensioned by (*) are not used in this algorithm
       dimension props(nprops), density(nblock),
     1 strainInc(nblock,ndir+nshr),
     2 relSpinInc(nblock,nshr), defgradOld(nblock,9),
     4 stressOld(nblock,ndir+nshr),
     5 stateOld(nblock,nstatev), enerInternOld(nblock),
     6 enerInelasOld(nblock),
     7 stretchNew(nblock,ndir+nshr), defgradNew(nblock,9),
     8 stressNew(nblock,ndir+nshr)
     
       dimension enerInelasNew(nblock),stateNew(nblock,nstatev),
     1 enerInternNew(nblock),eigVal(nblock,3)
     
       dimension nElement(nblock),nMatPoint(nblock),nLayer(nblock),
     1 nSecPoint(nblock)
     
      character*80 cmname
c      CHARACTER*8 CMNAME
C
c*****INPUT FORM*****************************************************************
!*USER MATERIAL,CONSTANTS=40
!#(-) Value in ax gives the exponent M for the yield function
!#       K         G        Cp        Xi        a1        a2        a3        a4
!        1         2         3         4         5         6         7         8
!#      a5        a6        a7        a8        a9       a10       a11       a12
!        9        10        11        12        13        14        15        16
!#Hardening
!        A         B         n    gt(JCX)   gs(JCX)   gc(JCX)   k(JCX)        (1) HARD PL
!       s0        Q1        C1        Q2        C2                            (2) HARD VOCE
!        A        e0         n        s0        Q1        C1     alpha        (3) HARD MSV
!       17        18        19        20        21        22        23        24
!#SRH+TS>JC
!# eps0dot         C         m        T0        Tr        Tm   epsAdot    0BE/1FE
!#      25        26        27        28        29        30        31        32
!#Failure criterion
! ADDFail (+) Tc, (-) EQPSmax
!      Wcr                           m12       m22       m44   ADDFail      (0/1) FAILURE CL(EQPS)
!        a         b         c         n        D4        D5   ADDFail      (2/3) FAILURE HCeps
!       D1        D2        D3        D4        D5             ADDFail      (4/5) FAILURE JCX
!        a         b         c       m12       m22       m44   ADDFail      (6/7) FAILURE HCsig
!      Wcr                           m12       m22       m44   ADDFail      (8/9) FAILURE CL(EQS)
!       33        34        35        36        37        38        39        40
!*DEPVAR, DELETE=7
!19
!1,EQPS,"Equivalent Plastic Strain"
!2,Seq,"Equivalent stress"
!3,Svm,"Equivalent VM stress"
!4,TRIAX,"Triaxiality"
!5,LODE,"Lode parameter"
!6,D,"Damage"
!7,FAIL,"Failure switch"
!8,EQS,"Equivalent strain"
!9,EQPSDot,"Equivalent Plastic Strain Rate"
!10,T,"Temperature"
!11,CHECK,"CHECK"
!12,ySRH,"Strain rate hardening"
!13,yTS,"Thermal softening"
!14,fSR,"Failure strain rate"
!15,fTS,"Failure thermal softening"
!16,Wcl,"CL plastic work"
!17,EQPSf,"Equivalent Plastic Strain to Failure"
!18,ITER,"Iterations"
!19,RES,"Residual"

      PARAMETER(ZERO=0.d0, ONE=1.d0, TWO=2.d0, THREE=3.d0,
     1  EIGHT=8.d0, FOUR=4.d0,
     2  THIRD=1.d0/3.d0, HALF=0.5d0, Op5=1.5d0, Od9=1.d0/9.d0)
C
      integer FAILFLAG, HARDFLAG
      INTEGER i, j, K1, K2, ii, NEWTON, Forflag, ITER
      real*8 PI, LAMBDA, KK, GG, G2, G3
      REAL*8 PPP(6,6), GGG(6,6), CCC(6,6), MM(6,6), L1(5,6), L2(5,6)
      REAL*8 S(6), St(6), Sd(6), N(6),SS1(6),SS2(6)
      REAL*8 CCxN(6)
      REAL*8 L1xS1(5), L2xS2(5), Sstar(6)
      REAL*8 NxCxN, check
      REAL*8 a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, my
      REAL*8 Sh, Seq, Svm, TRIAX, LODE
      REAL*8 eqpsOld, eqpsNew, dlambdaNew, dlambdaOld, ddlambda, deqps
      REAL*8 PHI, PHI1, PHI2, Mass, density, AepsV
      REAL*8 ySH, HSH, ySRH, HSRH, yTS, HTS, yOld, yNew
      REAL*8 epsVOld, epsVNew, PHIres
      REAL*8 AWp, AT, Xi, Cp, T0, Tr, Tm
      REAL*8 Dold, Dnew, FAILURE, fSR, ADDFAIL
      REAL*8 omega
      REAL*8 SIGMA1, SIGMA2, SIGMA3, PMAXEIGEN
      REAL*8 deqs, eqsOld, eqsNew, m12, m22, m44
      REAL*8 Eps0, EpsA, ySRH0, yTS0, fSR0, fTS0
C
C**********************************************************************
C     CONSTANTS
C**********************************************************************
      PI      = dacos(-1.d0)
      KK      = PROPS(1)
      GG      = PROPS(2)
      Cp      = PROPS(3)
      Xi      = PROPS(4)
      a1      = PROPS(5)
      a2      = PROPS(6)
      a3      = PROPS(7)
      a4      = PROPS(8)
      a5      = PROPS(9)
      a6      = PROPS(10)
      a7      = PROPS(11)
      a8      = PROPS(12)
      a9      = PROPS(13)
      a10     = PROPS(14)
      a11     = PROPS(15)
      a12     = PROPS(16)
      Eps0    = PROPS(25)
      C       = PROPS(26)
      T0      = PROPS(28)
      EpsA    = PROPS(31)
      Forflag = PROPS(32)
      m12     = PROPS(36)
      m22     = PROPS(37)
      m44     = PROPS(38)
      TOLER   = 1.0d-4
      NEWTON  = 150
c     HARDENING TYPE: (1) JCX; (2) VOCE; (3) MSV
      HARDflag = PROPS(24)
c     ADDITIONAL FAILURE
c     (0) NO FAILURE; (+) Critical Temperature; (-) MAX EQPS
      ADDFAIL  = PROPS(39)
c     FAILURE TYPE: (x/x) NO/YES
      FAILflag = PROPS(40)

c     Yld2000-3d m exponent
      IF (a1.lt.ZERO) THEN
        my = 1.d0
      END IF
      IF (a2.lt.ZERO) THEN
        my = 2.d0
      END IF
      IF (a3.lt.ZERO) THEN
        my = 3.d0
      END IF
      IF (a4.lt.ZERO) THEN
        my = 4.d0
      END IF
      IF (a5.lt.ZERO) THEN
        my = 5.d0
      END IF
      IF (a6.lt.ZERO) THEN
        my = 6.d0
      END IF
      IF (a7.lt.ZERO) THEN
        my = 7.d0
      END IF
      IF (a8.lt.ZERO) THEN
        my = 8.d0
      END IF
      IF (a9.lt.ZERO) THEN
        my = 9.d0
      END IF
      IF (a10.lt.ZERO) THEN
        my = 10.d0
      END IF
      IF (a11.lt.ZERO) THEN
        my = 11.d0
      END IF
      IF (a12.lt.ZERO) THEN
        my = 12.d0
      END IF

      a1    = abs(PROPS(5))
      a2    = abs(PROPS(6))
      a3    = abs(PROPS(7))
      a4    = abs(PROPS(8))
      a5    = abs(PROPS(9))
      a6    = abs(PROPS(10))
      a7    = abs(PROPS(11))
      a8    = abs(PROPS(12))
      a9    = abs(PROPS(13))
      a10   = abs(PROPS(14))
      a11   = abs(PROPS(15))
      a12   = abs(PROPS(16))

C**********************************************************************
C     COMPUTE THE MATRICES
C**********************************************************************
c     COMPUTE ELASTIC MODULI TENSOR
      G2     = TWO*GG
      LAMBDA = KK-G2/THREE

      DO i = 1, 6
        DO j = 1, 6
          CCC(i,j) = ZERO
        END DO
      END DO
      DO K1 = 1, 3
        DO K2 = 1, 3
          CCC(K2, K1) = LAMBDA
        END DO
        CCC(K1, K1) = G2 + LAMBDA
      END DO
      DO K1 = 4, 6
        CCC(K1, K1) = GG
      END DO
c     COMPUTE ANISOTROPIC FAILURE TENSOR
      DO i = 1, 6
        DO j = 1, 6
          MM(i,j) = ZERO
        END DO
      END DO
      MM(1,1) = ONE
      MM(1,2) = m12
      MM(2,2) = m22
      MM(3,3) = ONE
      MM(4,4) = m44
      MM(5,5) = ONE
      MM(6,6) = ONE
C     L1 ANISOTROPIC TENSOR
      DO i = 1, 5
        DO j = 1, 6
          L1(i,j) = ZERO
        END DO
      END DO
      L1(1,1) = THIRD*(TWO*a1)
      L1(1,2) = THIRD*(-a1)
      L1(1,3) = THIRD*(-a1)
      L1(2,1) = THIRD*(-a2)
      L1(2,2) = THIRD*(TWO*a2)
      L1(2,3) = THIRD*(-a2)
      L1(3,4) = a7
      L1(4,5) = a9
      L1(5,6) = a10

C     L2 ANISOTROPIC TENSOR
      DO i = 1, 5
        DO j = 1, 6
          L2(i,j) = ZERO
        END DO
      END DO
      L2(1,1) = Od9*(-TWO*a3+TWO*a4+EIGHT*a5-TWO*a6)
      L2(1,2) = Od9*(-FOUR*a4+FOUR*a6+a3-FOUR*a5)
      L2(1,3) = Od9*(a3+TWO*a4-FOUR*a5-TWO*a6)
      L2(2,1) = Od9*(FOUR*a3-FOUR*a4-FOUR*a5+a6)
      L2(2,2) = Od9*(-TWO*a3+EIGHT*a4+TWO*a5-TWO*a6)
      L2(2,3) = Od9*(-TWO*a3-FOUR*a4+TWO*a5+a6)
      L2(3,4) = a8
      L2(4,5) = a11
      L2(5,6) = a12


C**********************************************************************
C     Variable initialisation
C**********************************************************************

      IF ( stepTime .eq. zero ) THEN
        do k = 1, nblock
          call STRAIN_RATE_HARD(ySRH0,HSRH,Eps0,PROPS)
          call THERMAL_SOFT(yTS0,HSRH,T0,PROPS)
          call STRAIN_RATE_FAIL(fSR0,Eps0,PROPS)
          call THERMAL_SOFT_FAIL(fTS0,T0,PROPS)
          stateOld(k,9)   = Eps0           !INITIAL strain rate
          stateOld(k,10)  = T0             !INITIAL temperature
          stateOld(k,12)  = ySRH0          !INITIAL strain rate hardening
          stateOld(k,13)  = yTS0           !INITIAL thermal softening
          stateOld(k,14)  = fSR0           !INITIAL failure strain rate term
          stateOld(k,15)  = fTS0           !INITIAL failure thermal softening
          Aepsv = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)
           DO i= 1, NDIR
           stressNew(k,i)=stressOld(k,i)+LAMBDA*Aepsv+G2*strainInc(k,i)
           END DO
           DO i = NDIR+1, NDIR+NSHR
           stressNew(k,i)=stressOld(k,i)+G2*strainInc(k,i)
           END DO
        end do

      ELSE

        do k = 1, nblock
C**********************************************************************
C         STATE VARIABLES OLD
C**********************************************************************
          NUMEL   = nElement(k)
          eqpsOld = stateOld(k,1)     !EQPS*n+1
          Dold    = stateOld(k,6)     !Dn
          FAILURE = stateOld(k,7)     !FAILURE SWITCH
          eqsOld  = stateOld(k,8)     !Equivalent strain
          eqpsDot = stateOld(k,9)     !EQPSrate
          Told    = stateOld(k,10)    !Temperature
          check   = stateOld(k,11)    !check
          fSR     = stateOld(k,14)    !failure strain rate term
          fTS     = stateOld(k,15)    !failure thermal softening
          WclOld  = stateOld(k,16)    !CL plastic energy
          EQPSf   = stateOld(k,17)    !EQPSf
          ITER    = stateOld(k,18)    !IT to convergence
          PHIres  = stateOld(k,19)    !RESIDUAL
C**********************************************************************
C         ELASTIC PREDICTOR
C**********************************************************************
C         COMPUTE EQUIVALENT STRAIN
          IF (NSHR.eq.1) THEN
          strainInc(k,5) = ZERO
          strainInc(k,6) = ZERO
          ENDIF
          deqs   = sqrt(TWO/THREE*(
     1                  strainInc(k,1)**TWO+
     2                  strainInc(k,2)**TWO+
     3                  strainInc(k,3)**TWO+
     4             TWO*(strainInc(k,4)**TWO+
     5                  strainInc(k,5)**TWO+
     6                  strainInc(k,6)**TWO)))
          eqsNew = eqsOld + deqs

c         COMPUTE VOLUMETRIC STRAINe INCREMENT (Aepsv)
          Aepsv = strainInc(k,1)+strainInc(k,2)+strainInc(k,3)

c         COMPUTE TRIAL STRESS TENSOR (S*n+1=Sn+2G*Aee+K*Aepsv)
          DO i = 1, NDIR
          st(i) = stressOld(k,i) + LAMBDA*Aepsv + G2*strainInc(k,i)
          END DO
          DO i = NDIR+1, NDIR+NSHR
          st(i) = stressOld(k,i) + G2*strainInc(k,i)
          END DO

c         COMPUTE TRIAL HYDROSTATIC STRESS (Sh*n+1=Shn+1=Sh=1/3tr(s))
          sh = THIRD*(st(1)+st(2)+st(3))

c         COMPUTE TRIAL DEVIATORIC STRESS TENSOR (Sd*n+1=S*n+1-Sh*I)
          DO i = 1, NDIR
          sd(i) = st(i) - sh
          END DO
          DO i = NDIR+1, NDIR+NSHR
          sd(i) = st(i)
          END DO
          IF (NSHR.eq.1) THEN
          st(5) = ZERO
          st(6) = ZERO
          sd(5) = ZERO
          sd(6) = ZERO
          ENDIF

c         COMPUTE TRIAX AND INVARIANTS (SVM*n+1,TRIAX*n+1,LODE*n+1)
          CALL INVARIANTS(TRIAX,Svm,LODE,sd,SH)

C         COMPUTE YIELD FUNCTION Yld2000-3d
          CALL YLD_2000_3D(Seq,st,L1,L2,my,PROPS)

c**********************************************************************
c         HARDENING PL                                                *
c**********************************************************************
        IF (HARDflag.eq.1) THEN
c       Compute trial strain hardening Y*n+1 = function(EQPS*n+1)
        call PL_HARD(ySH,HSH,eqpsOld,PROPS)
        ENDIF
c**********************************************************************
c         HARDENING VOCE                                              *
c**********************************************************************
        IF (HARDflag.eq.2) THEN
c       Compute trial strain hardening Y*n+1 = function(EQPS*n+1)
        call VOCE_HARD(ySH,HSH,eqpsOld,PROPS)
        ENDIF
c**********************************************************************
c         HARDENING MIXED SWIFT-VOCE                                  *
c**********************************************************************
        IF (HARDflag.eq.3) THEN
c       Compute trial strain hardening Y*n+1 = function(EQPS*n+1)
        call SWIFT_VOCE_HARD(ySH,HSH,eqpsOld,PROPS)
        ENDIF
c**********************************************************************

c       Compute strain rate hardening
        call STRAIN_RATE_HARD(ySRH,HSRH,eqpsDot,PROPS)
c       Compute thermal softening
        call THERMAL_SOFT(yTS,HTS,TOld,PROPS)
c       TOTAL hardening
        yOld = ySH*ySRH*yTS

C**********************************************************************
C         YIELDING
C**********************************************************************
C         COMPUTE YIELD FUNCTION (phi=Seq*n+1-Y*n+1)
          PHI = seq - yOld

c         PLASTIC LOADING CONDITION
          IF ( PHI .gt. ZERO ) THEN 

C**********************************************************************
C         PLASTIC CORRECTOR (Backward-EULER/Forward-EULER)
C**********************************************************************
c         Al = phi / ( nT*C*n + H )
C         COMPUTE flow vector n
          CALL FLOW_VECTOR(N,st,L1,L2,my,PROPS)
          IF (NSHR.eq.1) THEN
          N(5) = ZERO
          N(6) = ZERO
          ENDIF

c         Compute C:n
          CALL PMxV(CCxN,CCC,N)
c         Compute n:C:n
          CALL PVxV(NxCxN,N,CCxN)

c         BACKWARD-EULER ALGORITHM
c         FIRST GUESS--------------------------------------------------
c         dAl_0
          dlambdaOld = ZERO
          ITER       = ZERO
c         ---------------------------------------------------------------

          DO KK = 1, NEWTON
            ddlambda = ZERO
            deqps    = ZERO
            AT       = ZERO
c           Plastic multiplier increment dAl for HC&JCX
   20       ddlambda = PHI/(NxCxN+
     1                 HSH*ySRH*yTS+
     2                 HSRH*(ONE/timeinc)*ySH*yTS+
     3                 HTS*(omega*Xi*Seq/(Cp*density(k)))*ySH*ySRH)
c           Plastic multiplier increment Al
   30       dlambdaNew = dlambdaOld + ddlambda
c           Equivalent plastic strain increment AEQPS
            deqps = dlambdaNew
C           ***********************************************************
C           PLASTIC UPDATE
C           ***********************************************************
c           Update EQPSk+1 (EQPSk+1 = EQPSn + Aeqps)
            eqpsNew = eqpsOld + deqps

c           Update EQPSrk+1
            eqpsDot = dmax1(Eps0,deqps/timeinc)

c           UPDATE THE STRESS TENSORS (Sdk+1=Sd*k+1-Al*C:m,Sn+1=Sk+1+Sh*I )
            DO i = 1, NDIR
            s(i)  = st(i) - deqps*CCxN(i)
            stressNew(k,i) = s(i)
            sd(i) = s(i) - sH
            ENDDO
            DO i = NDIR+1, NDIR+NSHR
            s(i)  = st(i) - deqps*CCxN(i)
            stressNew(k,i) = s(i)
            sd(i) = s(i)
            ENDDO
            IF (NSHR.eq.1) THEN
            s(5)  = ZERO
            s(6)  = ZERO
            sd(5) = ZERO
            sd(6) = ZERO
            ENDIF

C           UPDATE Seqn+1
            CALL YLD_2000_3D(Seq,s,L1,L2,my,PROPS)
        
c           Update Temperature Tk+1
            CALL OMEGA_LINEAR(omega,eqpsDot,PROPS)
            AWp  = Seq*deqps
            AT   = omega*Xi*AWp/(Cp*density(k))
            Tnew = Told + AT
           
c           **********************************************************************
c                    HARDENING PL                                                *
c           **********************************************************************
            IF (HARDflag.eq.1) THEN
c           Compute strain hardening Yk+1 = function(EQPSk+1)
            call PL_HARD(ySH,HSH,eqpsNew,PROPS)
            ENDIF
c           **********************************************************************
c                    HARDENING VOCE                                              *
c           **********************************************************************
            IF (HARDflag.eq.2) THEN
c           Compute strain hardening Yk+1 = function(EQPSk+1)
            call VOCE_HARD(ySH,HSH,eqpsNew,PROPS)
            ENDIF
c           **********************************************************************
c                    HARDENING MIXED SWIFT-VOCE                                  *
c           **********************************************************************
            IF (HARDflag.eq.3) THEN
c           Compute strain hardening Yk+1 = function(EQPSk+1)
            call SWIFT_VOCE_HARD(ySH,HSH,eqpsNew,PROPS)
            ENDIF
c           **********************************************************************
         
c           Compute strain rate hardening
            call STRAIN_RATE_HARD(ySRH,HSRH,eqpsDot,PROPS)
c           Compute thermal softening
            call THERMAL_SOFT(yTS,HTS,TNew,PROPS)
c           TOTAL hardening
            yNew = ySH*ySRH*yTS

c           CONVERGENCE CHECK******************************************
            PHI = Seq - yNew
      
      
c           Backward-Euler
            IF (DABS(PHI).LT.TOLER*yOld) GOTO 10
c           Forward-Euler
            IF (Forflag.eq.1) GOTO 10

            ITER = ITER + 1


C          WRITE WARNING MESSAGE TO THE .MSG FILE-------------------------
C          IF (ITER.EQ.NEWTON) THEN
C          Forflag = 1
C
C       write(*,*)
C       write(*,*)
C       write(*,*) '***************************************************'
C       write(*,*) '*                                                 *'
C       write(*,*) '*               W A R N I N G !                   *'
C       write(*,*) '*                                                 *'
C       write(*,*) '***************************************************'
C       write(*,*) '*                                                 *'
C       write(*,*) '* PLASTICITY ALGORITHM DID NOT CONVERGE AFTER     *'
C       write(*,*) '*',ITER,'ITERATIONS                          *'
C       write(*,*) '*-----------Forward-Euler solution----------------*'
C       write(*,*) '*                                                 *'
C       write(*,*) '***************************************************'
C       write(*,*)'ITER',ITER
C       write(*,*)'SH',HSH*ySRH*yTS*bOld
C       write(*,*)'SRH',HSRH*(ONE/timeinc)*ySH*yTS*bOld
C       write(*,*)'TS',HTS*(Xi*Seq/(Cp*density(k)))*ySH*ySRH*bOld
C       write(*,*)'dt',timeinc
C       write(*,*)'----------------------------------------------------'
C       write(*,*)'ddlambda',ddlambda,'dlambdaOld',dlambdaOld,
C     1            'dlambdaNew',dlambdaNew,'deqps',deqps
C       write(*,*)'eqpsNew',eqpsNew,'eqpsDot',eqpsDot,'TNew',TNew,
C     1            'DNew',DNew
C       write(*,*)'ySH',ySH,'ySRH',ySRH,'yTS',yTS,'bNew',bNew
C       write(*,*)'HSH',HSH,'HSRH',HSRH,'HTS',HTS,'Hbeta',Hbeta
C       write(*,*)'PHI',PHI,'TOLER',TOLER*PHI0,'yNew',yNew
C       write(*,*)'****************************************************'
C       STOP
C            IF (FAILflag.gt.1) GOTO 20
C            IF (FAILflag.le.1) GOTO 21
C
C           ENDIF
c           --------------------------------------------------------------



            IF (PHI.LT.ZERO) THEN
              IF ((ddlambda.GT.dlambdaOld).and.
     1            (dlambdaOld.GT.ZERO)) THEN
                ddlambda = dlambdaOld/TWO
                ELSE
                ddlambda = ddlambda/10.d0
              ENDIF
            GOTO 30
            ENDIF

            dLambdaOld = dLambdaNew

        END DO

   10  CONTINUE

C      RESIDUAL CHECK
       PHIres = DABS(PHI)/yOld

c        UPDATE TRIAX AND INVARIANTS (SVMn+1,TRIAXn+1,LODEn+1)
         CALL INVARIANTS(TRIAX,Svm,LODE,sd,SH)
c        **********************************************************************
c                 FAILURE CL(EQPS)                                            *
c        **********************************************************************
         IF (FAILflag.le.1) THEN
c        Compute S*=M:S
         CALL PMxV(Sstar,MM,s)
!        EIGENVALUES of S*
         CALL EIGENVALUES(SIGMA1,SIGMA2,SIGMA3,Sstar)
         PMAXEIGEN = dmax1(ZERO,SIGMA1)
c        Compute trial DAMAGE PARAMETER Dn+1
         call D_CL(Dnew,WclNew,DOld,WclOld,PMAXEIGEN,deqps,PROPS)
         fSR   = ONE
         fTS   = ONE
         EQPSf = ZERO
         ENDIF
c        **********************************************************************
c                 FAILURE HCeps                                               *
c        **********************************************************************
         IF ((FAILflag.eq.2).or.(FAILflag.eq.3)) THEN
c        Compute EQPSfn+1
         CALL FC_HC(fL,TRIAX,LODE,PROPS)
c        Compute failure strain rate term
         CALL STRAIN_RATE_FAIL(fSR,eqpsDot,PROPS)
c        Compute failure thermal softening term
         call THERMAL_SOFT_FAIL(fTS,TNew,PROPS)
c        TOTAL EQPSf
         EQPSf = fL*fSR*fTS
c        UPDATE DAMAGE PARAMETER Dn+1
         CALL D_LINEAR(Dnew,Dold,EQPSf,deqps,PROPS)
         WclNew = WclOld + Seq*deqps
         ENDIF
c        **********************************************************************
c                 FAILURE JCX                                                 *
c        **********************************************************************
         IF ((FAILflag.eq.4).or.(FAILflag.eq.5)) THEN
c        Compute  EQPSfn+1
         call FC_JCX(fL,TRIAX,LODE,PROPS)
c        Compute failure strain rate term
         CALL STRAIN_RATE_FAIL(fSR,eqpsDot,PROPS)
c        Compute failure thermal softening term
         call THERMAL_SOFT_FAIL(fTS,TNew,PROPS)
c        TOTAL EQPSf
         EQPSf = fL*fSR*fTS
c        UPDATE DAMAGE PARAMETER Dn+1
         CALL D_LINEAR(Dnew,Dold,EQPSf,deqps,PROPS)
         WclNew = WclOld + Seq*deqps
         ENDIF
c        **********************************************************************
c                 FAILURE HCsig                                               *
c        **********************************************************************
         IF ((FAILflag.eq.6).or.(FAILflag.eq.7)) THEN
c        Compute S*=M:S
         CALL PMxV(Sstar,MM,s)
!        EIGENVALUES of S*
         CALL EIGENVALUES(SIGMA1,SIGMA2,SIGMA3,Sstar)
c        Compute EQPSfn+1
         CALL FC_HCsig(Dnew,SIGMA1,SIGMA2,SIGMA3,PROPS)
         fSR   = ONE
         fTS   = ONE
         EQPSf = ZERO
         WclNew = WclOld + Seq*deqps
         ENDIF
c        **********************************************************************
c                 FAILURE CL(EQS)                                             *
c        **********************************************************************
         IF ((FAILflag.eq.8).or.(FAILflag.eq.9))  THEN
c        Compute S*=M:S
         CALL PMxV(Sstar,MM,s)
!        EIGENVALUES of S*
         CALL EIGENVALUES(SIGMA1,SIGMA2,SIGMA3,Sstar)
c        Compute trial DAMAGE PARAMETER Dn+1
         PMAXEIGEN = dmax1(ZERO,SIGMA1)
         call D_CL(Dnew,WclNew,DOld,WclOld,PMAXEIGEN,deqs,PROPS)
         fSR   = ONE
         fTS   = ONE
         EQPSf = ZERO
         ENDIF
c        **********************************************************************
         
          IF ((FAILFLAG.eq.0).or.
     1        (FAILFLAG.eq.2).or.
     2        (FAILFLAG.eq.4).or.
     3        (FAILFLAG.eq.6).or.
     4        (FAILFLAG.eq.8)) THEN
c           Elements when Damage = or > Dc
            IF (Dnew .GE. ONE) THEN
              Dnew = ONE
            END IF
          ENDIF

!=========FAILURE CRITERION ELEMENT REMOVAL============================
          IF ((FAILFLAG.eq.1).or.
     1        (FAILFLAG.eq.3).or.
     2        (FAILFLAG.eq.5).or.
     3        (FAILFLAG.eq.7).or.
     4        (FAILFLAG.eq.9)) THEN
c         Delete elements when Damage = or > Dc
          IF (Dnew .GE. ONE) THEN
            FAILURE = ZERO
            Dnew = ONE
          END IF
          ENDIF

!=========ADDITIONAL ELEMENT REMOVAL===================================
c         Delete elements when Temperature = or > Tc
          if (ADDFAIL.gt.ZERO)then
            if (Tnew.ge.ADDFAIL)then
            FAILURE = ZERO
            endif
          endif

c         Delete elements when EqPlasticStrain = or > MAXeqps
          if (ADDFAIL.lt.ZERO)then
            if (eqpsNew.gt.abs(ADDFAIL))then
            FAILURE = ZERO
            endif
          endif
!========================================================================

          enerInelasNew(k) = enerInelasOld(k) + AWp/density(k)
          check = ZERO

C         Update the state variables
          stateNew(k,1)  = eqpsNew        !EQPS
          stateNew(k,2)  = seq            !EQUIVALENT STRESS
          stateNew(k,3)  = sVM            !EQUIVALENT VM STRESS
          stateNew(k,4)  = TRIAX          !STRESS TRIAXIALITY
          stateNew(k,5)  = LODE           !LODE PARAMETER
          stateNew(k,6)  = Dnew           !DAMAGE PARAMETER
          stateNew(k,7)  = FAILURE        !FAILURE SWITCH
          stateNew(k,8)  = eqsNew         !Equivalent strain
          stateNew(k,9)  = eqpsDot        !Equivalent plastic strain rate
          stateNew(k,10) = Tnew           !Temperature
          stateNew(k,11) = omega          !omega
          stateNew(k,12) = ySRH           !Strain rate hardening
          stateNew(k,13) = yTS            !Thermal softening
          stateNew(k,14) = fSR            !Failure Strain rate
          stateNew(k,15) = fTS            !Failure thermal softening
          stateNew(k,16) = WclNew         !Cl plastic energy
          stateNew(k,17) = eqpsF          !EQPSf
          stateNew(k,18) = ITER           !IT to convergence
          stateNew(k,19) = PHIres         !RESIDUAL
		  
         IF (stateNew(k,1).GT.4.) THEN
          CALL XPLB_EXIT
         END IF

         ELSE
C**********************************************************************
C         ELASTIC UPDATE
C**********************************************************************
c        UPDATE THE STRESS TENSOR (Sdn+1=Sd*n+1,Sn+1=Sdn+1+Sh*I)
         DO i = 1, NDIR+NSHR
         s(i) = st(i)
         stressNew(k,i) = s(i)
         END DO
         IF (NSHR.eq.1) THEN
         s(5)  = ZERO
         s(6)  = ZERO
         sd(5) = ZERO
         sd(6) = ZERO
         ENDIF
         check   = ZERO
         omega   = ZERO
         Dnew    = Dold
         WclNew  = WclOld

         IF (DOld.ge.ONE) THEN
         Dnew = Dold
         GOTO 99
         ENDIF
c        ****************************************************************
c                 FAILURE HCsig                                         *
c        ****************************************************************
         IF ((FAILflag.eq.6).or.(FAILflag.eq.7)) THEN
c        Compute S*=M:S
         CALL PMxV(Sstar,MM,s)
!        EIGENVALUES of S*
         CALL EIGENVALUES(SIGMA1,SIGMA2,SIGMA3,Sstar)
c        Compute trial failure surface (n+1)*
         CALL FC_HCsig(Dnew,SIGMA1,SIGMA2,SIGMA3,PROPS)
         fSR   = ONE
         fTS   = ONE
         EQPSf = ZERO
         ENDIF
c        ****************************************************************
c                 FAILURE CL(EQS)                                       *
c        ****************************************************************
         IF ((FAILflag.eq.8).or.(FAILflag.eq.9))  THEN
c        Compute S*=M:S
         CALL PMxV(Sstar,MM,s)
!        EIGENVALUES of S*
         CALL EIGENVALUES(SIGMA1,SIGMA2,SIGMA3,Sstar)
c        Compute trial DAMAGE PARAMETER Dn+1*
         PMAXEIGEN = dmax1(ZERO,SIGMA1)
         call D_CL(Dnew,WclNew,DOld,WclOld,PMAXEIGEN,deqs,PROPS)
         fSR   = ONE
         fTS   = ONE
         EQPSf = ZERO
         ENDIF
c        ***************************************************************
   99    CONTINUE

         IF ((FAILFLAG.eq.6).or.
     1       (FAILFLAG.eq.8)) THEN
c          Elements when Damage = or > Dc
           IF (Dnew .GE. ONE) THEN
             Dnew = ONE
           END IF
         END IF

         IF  ((FAILFLAG.eq.7).or.
     1        (FAILFLAG.eq.9)) THEN
c          Delete elements when Damage = or > Dc
           IF (Dnew .GE. ONE) THEN
             FAILURE = ZERO
           END IF
         END IF

C        Update the state variables
         stateNew(k,1)  = eqpsOld        !EQPS
         stateNew(k,2)  = Seq            !EQUIVALENT STRESS
         stateNew(k,3)  = Svm            !EQUIVALENT SV STRESS
         stateNew(k,4)  = TRIAX          !STRESS TRIAXIALITY
         stateNew(k,5)  = LODE           !LODE PARAMETER
         stateNew(k,6)  = Dnew           !DAMAGE PARAMETER
         stateNew(k,7)  = FAILURE        !FAILURE SWITCH
         stateNew(k,8)  = eqsNew         !Equivalent strain
         stateNew(k,9)  = eqpsDot        !Equivalent plastic strain rate
         stateNew(k,10) = TOld           !Temperature
         stateNew(k,11) = omega          !omega function
         stateNew(k,12) = ySRH           !Strain rate hardening
         stateNew(k,13) = yTS            !Thermal softening
         stateNew(k,14) = fSR            !Strain rate hardening
         stateNew(k,15) = fTS            !Thermal softening
         stateNew(k,16) = WclNew         !CL plastic work
         stateNew(k,17) = eqpsF          !EQPSf
         stateNew(k,18) = ZERO           !IT to convergence
         stateNew(k,19) = ZERO           !RESIDUAL


        END IF

c          if (isnan(bOld))then
c          if(NUMEL.eq.48238)then
c          write(*,*)'============================================'
c          write(*,*)'time',totaltime,TRIAX,LODE,EQPSF,Dnew,bNew
c          endif

        end do
      end if
C
      return
      end















CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     AUXILIARY SUBROUTINES                                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INVARIANTS(TRIAX,SEQ,LODE,VECTOR,SH)
c     Calculate invariants of a symmetric matrix
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 J2, J3, CHI, TRIAX, SEQ, LODE, SH, PI, VECTOR(6)
      PARAMETER(zero=0.d0,one=1.d0,two=2.d0,three = 3.d0,tseven= 27.d0,
     1          HALF = 0.5d0)
      PI = dacos(-ONE)

      J2 = HALF*
     1   (VECTOR(1)*VECTOR(1)+VECTOR(2)*VECTOR(2)+VECTOR(3)*VECTOR(3))
     1   +VECTOR(4)*VECTOR(4)+VECTOR(5)*VECTOR(5)+VECTOR(6)*VECTOR(6)
      J3 = VECTOR(1)*VECTOR(2)*VECTOR(3)
     1    +TWO*VECTOR(4)*VECTOR(5)*VECTOR(6)
     2    -VECTOR(1)*VECTOR(5)*VECTOR(5)
     3    -VECTOR(2)*VECTOR(6)*VECTOR(6)
     4    -VECTOR(3)*VECTOR(4)*VECTOR(4)

      IF (J2.eq.ZERO) THEN
        SEQ   = ZERO
        TRIAX = ZERO
        CHI   = ZERO
      ELSE
        SEQ   = DSQRT(THREE*J2)
        TRIAX = SH/SEQ
        CHI   = (TSEVEN*J3)/(TWO*SEQ**THREE)
        IF (CHI.gt.ONE ) THEN
          CHI = ONE
        END IF
        IF(CHI.lt.-ONE) THEN
          CHI = -ONE
        END IF
      END IF
      LODE   = ONE - (TWO/PI)*DACOS(CHI)
      RETURN
      END

      SUBROUTINE PMXV(VECOUT,MAT,VEC)
C     Matrix product
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER i,j
      DIMENSION MAT(6,6),VEC(6),VECOUT(6)
      DO i = 1, 6
        VECOUT(i) = 0.D0
        DO j= 1, 6
        VECOUT(i) = VECOUT(i) + MAT(i,j)*VEC(j)
        END DO
      END DO
      RETURN
      END

      SUBROUTINE PMXV56(VECOUT,MAT,VEC)
C     Matrix product
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER ii, jj
      DIMENSION MAT(5,6),VEC(6),VECOUT(5)
      DO ii = 1, 5
        VECOUT(ii) = 0.D0
        DO jj= 1, 6
        VECOUT(ii) = VECOUT(ii) + MAT(ii,jj)*VEC(jj)
        END DO
      END DO
      RETURN
      END

      SUBROUTINE PVxV(POUT,VEC1,VEC2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION VEC1(6),VEC2(6)
      REAL*8 POUT
      INTEGER i,j
      POUT = 0.d0
        DO i=1,6
        POUT = POUT + VEC1(i)*VEC2(i)
        END DO
      RETURN 
      END

      SUBROUTINE PCVxV(POUT,VEC1,VEC2)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION VEC1(6),VEC2(6)
      REAL*8 POUT
      INTEGER i,j
      POUT = 0.d0
      TWO  = 2.D0
        DO i=1,3
        POUT = POUT + VEC1(i)*VEC2(i)
        END DO
        DO i=4,6
        POUT = POUT + TWO*VEC1(i)*VEC2(i)
        END DO
      RETURN 
      END

      SUBROUTINE SWIFT_VOCE_HARD(YSH,HSH,eqps,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 A, e0, n, Q, b, s0, alpha, HSH, YSH, eqps
      ONE=1.D0
      ! k=alpha*k1(eps,p)+(1-alpha)*k2(eps,p)
      ! k1=A*(eps,p+e0)**n
      ! k2=Q*(1-exp(-b*eps,p))+sig0
      A     = PROPS(17)
      e0    = PROPS(18)
      n     = PROPS(19)
      s0    = PROPS(20)
      Q     = PROPS(21)
      C     = PROPS(22)
      alpha = PROPS(23)
      YSH=alpha*A*(eqps+e0)**n+(ONE-alpha)*(s0+Q*(ONE-EXP(-C*eqps)))
      HSH=alpha*A*n*(eqps+e0)**(n-ONE)+(ONE-alpha)*Q*C*EXP(-C*eqps)
      RETURN
      END 

      SUBROUTINE VOCE_HARD(YSH,HSH,eqps,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 s0, Q1, C1, Q2, C2, HSH, YSH, eqps
      ONE=1.D0
      ! Y=s0+Q1*(1-exp(-C1*eps))+Q2*(1-exp(-C2*eps))
      s0    = PROPS(17)
      Q1    = PROPS(18)
      C1    = PROPS(19)
      Q2    = PROPS(20)
      C2    = PROPS(21)
      YSH=s0+Q1*(ONE-EXP(-C1*eqps))+Q2*(ONE-EXP(-C2*eqps))
      HSH=Q1*C1*EXP(-C1*eqps)+Q2*C2*EXP(-C2*eqps)
      RETURN
      END 

      SUBROUTINE PL_HARD(YSH,HSH,eqps,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 s0, B, n, HSH, YSH, eqps
      ZERO=0.D0
      ONE=1.D0
      ! Y = s0+B*eps^n
      s0    = PROPS(17)
      B     = PROPS(18)
      n     = PROPS(19)
      IF (eqps.le.ZERO) THEN
        YSH = s0
        HSH = ZERO
      ELSE
        YSH = s0+B*eqps**n
        HSH = n*B*eqps**(n-ONE)
      ENDIF
      RETURN
      END

      SUBROUTINE STRAIN_RATE_HARD(ySRH,HSRH,eqpsDot,PROPS)
      include 'vaba_param.inc'
      DIMENSION PROPS(*)
      REAL*8 ySRH, C, eqpsDot, eps0
      ZERO = 0.D0
      ONE  = 1.D0
      ! ySRH  = 1+C*ln(eqpsDot/eps0)
      C     = PROPS(26)
      Eps0  = PROPS(25)
      IF ((eqpsDot .le. Eps0).or.(C.eq.ZERO)) THEN
        ySRH = ONE
        HSRH = ZERO
        ELSE
        ySRH = ONE + C*log(eqpsDot/Eps0)
        HSRH = C/eqpsDot
      ENDIF
      RETURN
      END

      SUBROUTINE THERMAL_SOFT(yTS,HTS,Temp,PROPS)
      include 'vaba_param.inc'
      DIMENSION PROPS(*)
      REAL*8 yTS, Temp, T0, Tr, Tm, m, Tstar
      ZERO = 0.D0
      ONE  = 1.D0
      m      = PROPS(27)
      T0     = PROPS(28)
      Tr     = PROPS(29)
      Tm     = PROPS(30)
      ! yTS  = (1-(T*)^m)
      ! T*   = (T-Tr)/(Tm-Tr)
      IF ((Temp .le. Tr).or.(m.eq.ZERO)) THEN
        yTS = ONE
        HTS = ZERO
        ELSE
        Tstar = (Temp-Tr)/(Tm-Tr)
        yTS = ONE - Tstar**m
        HTS = -ONE/(Tm-Tr)*m*Tstar**(m-ONE)
      ENDIF
      RETURN
      END

      SUBROUTINE OMEGA_LINEAR(omega,eqpsDot,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8  omega, eqpsDot
      REAL*8  denom, Eps0, EpsA
      REAL*8  a_par(4)
      PARAMETER (ZERO = 0.D0, ONE  = 1.D0, TWO = 2.D0, THREE = 3.D0)
      Eps0    = PROPS(25)
      EpsA    = PROPS(31)
      IF (eqpsDot .LE. Eps0) THEN
        omega  = ZERO
        ELSE IF (eqpsDot .GT. EpsA) THEN
        omega  = ONE  
        ELSE 
        denom     = ( Eps0-EpsA )**THREE
        a_par(0)  = ( Eps0**TWO * (-THREE*EpsA+Eps0) ) / denom
        a_par(1)  = ( TWO*THREE*EpsA*Eps0 ) / denom
        a_par(2)  = -( THREE*(Eps0+EpsA) ) / denom
        a_par(3)  =   TWO / denom
        omega = a_par(0)+a_par(1)*eqpsDot+
     1        a_par(2)*(eqpsDot)**TWO+a_par(3)*(eqpsDot)**THREE
      END IF
      RETURN
      END

      SUBROUTINE FC_HC(EQPSf,TRIAX,LODE,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 EQPSf, TRIAX, LODE
      REAL*8 a, b, c, n, f1, f2, f3, PI
      PARAMETER (ZERO = 0.D0, HALF = 0.5D0, ONE  = 1.D0, OP5 = 1.5D0, 
     1           TWO = 2.D0, THREE = 3.D0, SIX = 6.D0)
      PI = DACOS(-ONE)
      a = PROPS(33)
      b = PROPS(34)
      c = PROPS(35)
      n = PROPS(36)
      f1 =  (TWO/THREE)*DCOS((PI/SIX)*(ONE-LODE))
      f2 =  (TWO/THREE)*DCOS((PI/SIX)*(THREE+LODE))
      f3 = -(TWO/THREE)*DCOS((PI/SIX)*(ONE+LODE))
      EQPSf = (b*(ONE+c)**(ONE/n))*
     1        ((((HALF*
     2        ((dabs(f1-f2)**a)+(dabs(f2-f3)**a)+
     2         (dabs(f1-f3)**a)))**(ONE/a))+
     3        c*(TWO*TRIAX+f1+f3))**(-ONE/n))
      RETURN
      END

      SUBROUTINE FC_HCsig(D,SIGMA1,SIGMA2,SIGMA3,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8  D, SIGMA1, SIGMA2, SIGMA3
      REAL*8 Shc, a, b, c
      PARAMETER (ZERO = 0.D0, HALF = 0.5D0, ONE  = 1.D0)
      a = PROPS(33)
      b = PROPS(34)
      c = PROPS(35)
      Shc = (HALF*(
     1             (SIGMA1-SIGMA2)**a+(SIGMA2-SIGMA3)**a+
     2             (SIGMA1-SIGMA3)**a))**(ONE/a)
      
      D = (Shc + c*(SIGMA1+SIGMA3))/b
      RETURN
      END

      SUBROUTINE FC_JCX(EQPSf,TRIAX,LODE,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 EQPSf, TRIAX, LODE
      REAL*8 D1, D2, D3, D4, D5, gt, gs, gc, MU, Lk
      gt = PROPS(20)
      gs = PROPS(21)
      gc = PROPS(22)
      Lk = PROPS(23)
      D1 = PROPS(33)
      D2 = PROPS(34)
      D3 = PROPS(35)
      IF (LODE.gt.ZERO) THEN
        MU = gs + (gt-gs)*(dabs(LODE))**Lk
        ELSE
        MU = gs + (gc-gs)*(dabs(LODE))**Lk
      ENDIF
      EQPSf = (D1+D2*dEXP(D3*TRIAX))*MU
      RETURN
      END

      SUBROUTINE D_CL(Dnew,WclNew,Dold,WclOld,SIGMA1,deqps,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 SIGMA1, deqps
      REAL*8 Dnew, Dold, WCR, WclNew, WclOld, dWcl, dD
      WCR    = PROPS(33)
      dWcl   = SIGMA1*deqps
      dD     = dWcl/WCR
      WclNew = WclOld + dWcl
      Dnew   = Dold   + dD
      RETURN
      END

      SUBROUTINE STRAIN_RATE_FAIL(fSR,eqpsDot,PROPS)
      include 'vaba_param.inc'
      DIMENSION PROPS(*)
      REAL*8 fSR, D4, eqpsdot, eps0
      ZERO = 0.D0
      ONE  = 1.D0
      ! fSR  = 1+D4*ln(EQPSDot/eps0)
      Eps0  = PROPS(25)
      D4    = PROPS(37)
      IF ((D4.eq.ZERO).or.(eqpsdot .le. Eps0)) THEN
        fSR = ONE
        ELSE
        fSR = ONE + D4*log(eqpsDot/Eps0)
      ENDIF
      RETURN
      END

      SUBROUTINE THERMAL_SOFT_FAIL(fTS,Temp,PROPS)
      include 'vaba_param.inc'
      DIMENSION PROPS(*)
      REAL*8 fTS, Temp, T0, Tr, Tm, Tstar, D5
      ZERO = 0.D0
      ONE  = 1.D0
      T0     = PROPS(28)
      Tr     = PROPS(29)
      Tm     = PROPS(30)
      D5     = PROPS(38)
      ! fTS  =  1+D5(T*)
      ! T*   = (T-Tr)/(Tm-Tr)
      IF ((Temp .le. Tr).or.(D5.eq.ZERO)) THEN
        fTS = ONE
        ELSE
        Tstar = (Temp-Tr)/(Tm-Tr)
        fTS = ONE - D5*Tstar
      ENDIF
      RETURN
      END

      SUBROUTINE D_LINEAR(Dnew,Dold,EQPSf,deqps,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 Dold, Dnew, dD, deqps, ONE
      ZERO = 0.D0
      ONE  = 1.D0
      dD = dmax1(ZERO, (ONE/EQPSf)*deqps)
      Dnew = Dold + dD
      RETURN
      END

!     3d extension of the YLD2000-2d
      SUBROUTINE YLD_2000_3D(Seq,sigma,L1,L2,my,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      INTEGER ii
      REAL*8 sigma(6),L1(5,6),L2(5,6),L1xS(5),L2xS(5),SS1(6),SS2(6)
      REAL*8 PHI1, PHI2, Seq, my
      PARAMETER( zero = 0.d0, one = 1.d0, two = 2.d0, THREE = 3.d0,
     1  FOUR = 4.d0, half = 0.5d0)

      CALL PMXV56(L1xS,L1,sigma)
      CALL PMXV56(L2xS,L2,sigma)

      SS1(1) = L1xS(1)
      SS1(2) = L1xS(2)
      SS1(3) = - L1xS(1) - L1xS(2)
      SS1(4) = L1xS(3)
      SS1(5) = L1xS(4)
      SS1(6) = L1xS(5)
      SS2(1) = L2xS(1)
      SS2(2) = L2xS(2)
      SS2(3) = - L2xS(1) - L2xS(2)
      SS2(4) = L2xS(3)
      SS2(5) = L2xS(4)
      SS2(6) = L2xS(5)
          
      PHI1 = ((SS1(1)-SS1(2))**TWO+
     1        FOUR*(SS1(4)**TWO+SS1(5)**TWO+SS1(6)**TWO))**(HALF*my)
      PHI2 = (THREE/TWO*(SS2(1)+SS2(2))+
     1        HALF*dsqrt((SS2(1)-SS2(2))**TWO+
     2        FOUR*(SS2(4)**TWO+SS2(5)**TWO+SS2(6)**TWO)))**my+
     3       (THREE/TWO*(SS2(1)+SS2(2))-
     4        HALF*dsqrt((SS2(1)-SS2(2))**TWO+
     5        FOUR*(SS2(4)**TWO+SS2(5)**TWO+SS2(6)**TWO)))**my
      Seq  = (HALF*(PHI1+PHI2))**(ONE/my)
      RETURN
      END


      SUBROUTINE FLOW_VECTOR(N,stress,L1,L2,my,PROPS)
      IMPLICIT DOUBLE PRECISION (A-Z)
      DIMENSION PROPS(*)
      REAL*8 L1(5,6),L2(5,6)
      REAL*8 stress(6),stressDs(6),N(6)
      REAL*8 Seq, my, SeqOld, SeqNew, ds
      ds = 5.0D-6

C     Equivalent stress
      CALL YLD_2000_3D(Seq,stress,L1,L2,my,PROPS)
      SeqOld = Seq
     
C     N1 Equivalent stress+Dstress1
      stressDs(1)=stress(1)+ds
      stressDs(2)=stress(2)
      stressDs(3)=stress(3)
      stressDs(4)=stress(4)
      stressDs(5)=stress(5)
      stressDs(6)=stress(6)
      CALL YLD_2000_3D(Seq,stressDs,L1,L2,my,PROPS)
      SeqNew = Seq
      N(1)= (SeqNew-SeqOld)/ds
C     
C     N2 Equivalent stress+Dstress2
      stressDs(1)=stress(1)
      stressDs(2)=stress(2)+ds
      stressDs(3)=stress(3)
      stressDs(4)=stress(4)
      stressDs(5)=stress(5)
      stressDs(6)=stress(6)
      CALL YLD_2000_3D(Seq,stressDs,L1,L2,my,PROPS)
      SeqNew = Seq
      N(2)= (SeqNew-SeqOld)/ds

C     N3 Equivalent stress+Dstress3
      stressDs(1)=stress(1)
      stressDs(2)=stress(2)
      stressDs(3)=stress(3)+ds
      stressDs(4)=stress(4)
      stressDs(5)=stress(5)
      stressDs(6)=stress(6)
      CALL YLD_2000_3D(Seq,stressDs,L1,L2,my,PROPS)
      SeqNew = Seq
      N(3)= (SeqNew-SeqOld)/ds

C     N4 Equivalent stress+Dstress4
      stressDs(1)=stress(1)
      stressDs(2)=stress(2)
      stressDs(3)=stress(3)
      stressDs(4)=stress(4)+ds
      stressDs(5)=stress(5)
      stressDs(6)=stress(6)
      CALL YLD_2000_3D(Seq,stressDs,L1,L2,my,PROPS)
      SeqNew = Seq
      N(4)= (SeqNew-SeqOld)/ds
C     
C     N5 Equivalent stress+Dstress5
      stressDs(1)=stress(1)
      stressDs(2)=stress(2)
      stressDs(3)=stress(3)
      stressDs(4)=stress(4)
      stressDs(5)=stress(5)+ds
      stressDs(6)=stress(6)
      CALL YLD_2000_3D(Seq,stressDs,L1,L2,my,PROPS)
      SeqNew = Seq
      N(5)= (SeqNew-SeqOld)/ds
C     
C     N6 Equivalent stress+Dstress6
      stressDs(1)=stress(1)
      stressDs(2)=stress(2)
      stressDs(3)=stress(3)
      stressDs(4)=stress(4)
      stressDs(5)=stress(5)
      stressDs(6)=stress(6)+ds
      CALL YLD_2000_3D(Seq,stressDs,L1,L2,my,PROPS)
      SeqNew = Seq
      N(6)= (SeqNew-SeqOld)/ds

      CONTINUE
      RETURN
      END

      SUBROUTINE EIGENVALUES(SIGMA1,SIGMA2,SIGMA3,S)
      IMPLICIT DOUBLE PRECISION (A-Z)
      REAL*8 SIGMA1, SIGMA2, SIGMA3
      REAL*8 S(6), Sdev(6)
      REAL*8 Svm, LODE, TRIAX, f1, f2, f3
      PARAMETER (ZERO = 0.D0, ONE  = 1.D0, 
     1           TWO = 2.D0, THREE = 3.D0, SIX = 6.D0)
      PI = dacos(-ONE)
c     GET THE HYDROSTATIC STRESS AND DEVIATOR
      SH = ONE/THREE*(S(1)+S(2)+S(3))
      DO j = 1,3
      Sdev(j) = S(j) - SH
      ENDDO
      DO j = 4,6
      Sdev(j) = S(j)
      ENDDO
c     GET THE INVARIANTS
      CALL INVARIANTS(TRIAX,Svm,LODE,Sdev,SH)
      f1 =  (TWO/THREE)*DCOS((PI/SIX)*(ONE-LODE))
      f2 =  (TWO/THREE)*DCOS((PI/SIX)*(THREE+LODE))
      f3 = -(TWO/THREE)*DCOS((PI/SIX)*(ONE+LODE))
c     Compute eigenvalues of TENSOR
      SIGMA1 = Svm*(TRIAX+f1)
      SIGMA2 = Svm*(TRIAX+f2)
      SIGMA3 = Svm*(TRIAX+f3)
      RETURN
      END
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


