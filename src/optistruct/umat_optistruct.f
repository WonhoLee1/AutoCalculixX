C ============================================================
C  WHTOOLs MATCALIB 2026 — OptiStruct Material User Subroutine
C  PRF (Parallel Rheological Framework) for solids (C3D8 etc.)
C
C  Build (Intel Fortran on Windows):
C    ifort /dll /libs:dll /threads umat_optistruct.f /Fe:umat.dll
C
C  Build (gfortran on Windows with MSYS2):
C    gfortran -shared -o umat.dll umat_optistruct.f
C
C  Build (gfortran on Linux):
C    gfortran -fPIC -shared -o umat.so umat_optistruct.f
C
C  OptiStruct input (MATUSR):
C    MATUSR, MID, USUBID, NDEPVAR, GROUP, DENSITY, PROPERTY
C            PROPi = N_NETWORK, FLAG_HE, FLAG_PL, C10, C01, ...
C
C  USUBID: user subroutine ID to distinguish material types
C          (passed as idu parameter)
C ============================================================
C
C  REQUIRED subroutines:
C    1. usermaterial  — nonlinear large displacement
C    2. smatusr       — linear / small displacement (stub)
C    3. initusr       — state variable label initialization (optional)
C
C  INCLUDED subroutines (from CalculiX UMAT):
C    PRF_EVAL, PRF_UPDATE, PRF_INIT_STATE, HYPER_STRESS
C    POLY_STRESS, ARRUDA_BOYCE, YEOH_STRESS, GENT_STRESS,
C    OGDEN_STRESS, JACOBI3
C    VISC_BB, VISC_SINH, VISC_POWER
C    KMATINV3, PRODMAT, PRODAAT, MATMUL_3X3, INV3X3, CALC_MATB
C ============================================================


C ============================================================
C  usermaterial — Nonlinear large displacement analysis
C  Called at each integration point, each iteration
C ============================================================
      SUBROUTINE usermaterial(
     &     idu,                                      ! user material ID
     &     stress,                                   ! stress (in/out)
     &     dstrain,                                  ! strain increment
     &     dfgrd0, dfgrd1,                           ! F_old, F_new
     &     statev,                                   ! state variables (in/out)
     &     nstatev,                                  ! number of state vars
     &     drot,                                     ! rotation increment
     &     ddsdde,                                   ! tangent stiffness (out)
     &     ntens, ndi, nshear,                       ! tensor dimensions
     &     nprops, props,                            ! material properties
     &     strain,                                   ! total strain
     &     time, dtime,                              ! time info
     &     temp, dtemp,                              ! temperature
     &     predef, dpred,                            ! predefined fields
     &     coords,                                   ! coordinates
     &     celent,                                   ! characteristic length
     &     noel, npt, layer, kspt,                   ! element info
     &     kinc, kstep,                              ! increment/step
     &     ieuid,                                    ! element ID
     &     ierr)                                     ! error flag
C
      IMPLICIT NONE
C
C ---- arguments -------------------------------------------------
      INTEGER idu, nstatev, ntens, ndi, nshear, nprops
      INTEGER noel, npt, layer, kspt, kinc, kstep, ieuid
      DOUBLE PRECISION stress(ntens), dstrain(ntens)
      DOUBLE PRECISION dfgrd0(3,3), dfgrd1(3,3)
      DOUBLE PRECISION statev(nstatev), drot(3,3)
      DOUBLE PRECISION ddsdde(ntens,ntens), props(nprops)
      DOUBLE PRECISION strain(ntens), time(2), dtime
      DOUBLE PRECISION temp, dtemp, predef(*), dpred(*)
      DOUBLE PRECISION coords(3), celent
      INTEGER ierr
C
C ---- local variables -------------------------------------------
      DOUBLE PRECISION F(3,3), SIG_CAUCHY(3,3), PK2(3,3)
      DOUBLE PRECISION SIG_PERT(6), SIG_REF(6)
      DOUBLE PRECISION FINV(3,3), FINV_T(3,3), VEC(3,3)
      DOUBLE PRECISION FPERT(3,3), DELTAF(3,3)
      DOUBLE PRECISION F_SAVE(3,3)
      DOUBLE PRECISION AJ, DET, EPS, PERT
      INTEGER NTERMS, N_NETWORK, FLAG_HE, FLAG_PL, NSTATV
      INTEGER NHE, I, J, K, ICOL, P, Q
      INTEGER NPOLY, ISTART, IP
      PARAMETER (EPS=1.0D-8)
C
C     Voigt index mapping: 1→11, 2→22, 3→33, 4→12, 5→13, 6→23
      INTEGER VPMAP(2,6)
      DATA VPMAP /1,1, 2,2, 3,3, 1,2, 1,3, 2,3/
C
      DOUBLE PRECISION ZERO, ONE, TWO, THREE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
C
C ---- parse material properties ---------------------------------
C     PROPS(1) = N_NETWORK
C     PROPS(2) = FLAG_HE
C     PROPS(3) = FLAG_PL
C     PROPS(4..) = hyperelastic + viscous + G, RBULK
C
      N_NETWORK = NINT(PROPS(1))
      FLAG_HE   = NINT(PROPS(2))
      FLAG_PL   = NINT(PROPS(3))
      NTERMS    = nprops
      NSTATV    = 2 + 12 * N_NETWORK
C
      IF (NSTATV .GT. nstatev) THEN
        IERR = 1
        RETURN
      END IF
C
C ---- copy deformation gradient (OptiStruct: dfgrd = F) --------
      DO I = 1, 3
        DO J = 1, 3
          F(I,J) = dfgrd1(I,J)
        END DO
      END DO
C
C ---- initialize state variables (first call) -------------------
      IF (time(1) .LT. EPS .OR. STATEV(1) .LT. EPS) THEN
        CALL PRF_INIT_STATE_OS(STATEV, N_NETWORK, FLAG_HE, PROPS(4))
      END IF
C
C ---- 5. compute J = det(F) --------------------------------------
      AJ = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
     1   - F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
     2   + F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
      STATEV(2) = AJ
C
C ---- 6. update internal variables (creep update) ----------------
      CALL PRF_UPDATE(F, DTIME, PROPS, NTERMS, N_NETWORK,
     &                FLAG_HE, FLAG_PL, STATEV, NSTATV)
C
C ---- 7. compute total Cauchy stress -----------------------------
      CALL SIGMA(F, DTIME, PROPS, NTERMS, N_NETWORK, FLAG_HE,
     &           FLAG_PL, STATEV, NSTATV, SIG_CAUCHY)
C
C ---- 8. Cauchy to PK2: S = J * F^-1 * sigma * F^-T -------------
      CALL INV3X3(F, FINV, DET)
      DO I = 1, 3
        DO J = 1, 3
          FINV_T(I,J) = FINV(J,I)
        END DO
      END DO
      CALL MATMUL_3X3(FINV, SIG_CAUCHY, VEC)
      CALL MATMUL_3X3(VEC, FINV_T, PK2)
      DO I = 1, 3
        DO J = 1, 3
          PK2(I,J) = AJ * PK2(I,J)
        END DO
      END DO
C
C ---- store stress in Voigt notation -----------------------------
      STRESS(1) = SIG_CAUCHY(1,1)
      STRESS(2) = SIG_CAUCHY(2,2)
      STRESS(3) = SIG_CAUCHY(3,3)
      STRESS(4) = SIG_CAUCHY(1,2)
      IF (ntens .GE. 5) STRESS(5) = SIG_CAUCHY(1,3)
      IF (ntens .GE. 6) STRESS(6) = SIG_CAUCHY(2,3)
C
C ---- save reference PK2 for tangent -----------------------------
      SIG_REF(1) = PK2(1,1)
      SIG_REF(2) = PK2(2,2)
      SIG_REF(3) = PK2(3,3)
      SIG_REF(4) = PK2(1,2)
      SIG_REF(5) = PK2(1,3)
      SIG_REF(6) = PK2(2,3)
C
C ---- numerical tangent stiffness (Miehe perturbation) -----------
      DO I = 1, ntens
        DO J = 1, ntens
          DDSDDE(I,J) = ZERO
        END DO
      END DO
C
      DO I = 1, 3
        DO J = 1, 3
          FINV_T(I,J) = FINV(J,I)
        END DO
      END DO
C
      DO ICOL = 1, ntens
        P = VPMAP(1, ICOL)
        Q = VPMAP(2, ICOL)
C
C       perturbation: δF = eps * F^-T * (e_p ⊗ e_q)
        DO I = 1, 3
          DO J = 1, 3
            DELTAF(I,J) = ZERO
          END DO
        END DO
        IF (P .EQ. Q) THEN
          DO I = 1, 3
            DELTAF(I,P) = EPS * FINV_T(I,P)
          END DO
        ELSE
          DO I = 1, 3
            DELTAF(I,Q) = EPS * FINV_T(I,P)
          END DO
        END IF
C
C       save state before perturbation (for frozen tangent)
        DO K = 1, NSTATV
          F_SAVE(1,K) = STATEV(K)
        END DO
C
        DO I = 1, 3
          DO J = 1, 3
            FPERT(I,J) = F(I,J) + DELTAF(I,J)
          END DO
        END DO
C
C       restore state (frozen Fv tangent)
        DO K = 1, NSTATV
          STATEV(K) = F_SAVE(1,K)
        END DO
C
C       stress at perturbed F
        CALL SIGMA(FPERT, DTIME, PROPS, NTERMS, N_NETWORK,
     &             FLAG_HE, FLAG_PL, STATEV, NSTATV, SIG_CAUCHY)
C
C       Cauchy to PK2
        CALL INV3X3(FPERT, FINV, DET)
        DO I = 1, 3
          DO J = 1, 3
            FINV_T(I,J) = FINV(J,I)
          END DO
        END DO
        DET = FPERT(1,1)*(FPERT(2,2)*FPERT(3,3)
     1       -FPERT(2,3)*FPERT(3,2))
     2      - FPERT(1,2)*(FPERT(2,1)*FPERT(3,3)
     3       -FPERT(2,3)*FPERT(3,1))
     4      + FPERT(1,3)*(FPERT(2,1)*FPERT(3,2)
     5       -FPERT(2,2)*FPERT(3,1))
        CALL MATMUL_3X3(FINV, SIG_CAUCHY, VEC)
        CALL MATMUL_3X3(VEC, FINV_T, PK2)
        DO I = 1, 3
          DO J = 1, 3
            PK2(I,J) = DET * PK2(I,J)
          END DO
        END DO
C
        SIG_PERT(1) = PK2(1,1)
        SIG_PERT(2) = PK2(2,2)
        SIG_PERT(3) = PK2(3,3)
        SIG_PERT(4) = PK2(1,2)
        IF (ntens .GE. 5) SIG_PERT(5) = PK2(1,3)
        IF (ntens .GE. 6) SIG_PERT(6) = PK2(2,3)
C
C       finite difference
        DO I = 1, ntens
          DDSDDE(I,ICOL) = (SIG_PERT(I) - SIG_REF(I)) / EPS
        END DO
      END DO
C
      IERR = 0
      RETURN
      END


C ============================================================
C  SIGMA — compute Cauchy stress (wraps PRF_EVAL)
C ============================================================
      SUBROUTINE SIGMA(F, DT, ELCONLOC, NTERMS,
     &                 N_NETWORK, FLAG_HE, FLAG_PL,
     &                 STATEV, NSTATV, SIG)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), DT, ELCONLOC(*), STATEV(*), SIG(3,3)
      INTEGER NTERMS, N_NETWORK, FLAG_HE, FLAG_PL, NSTATV
C
      DOUBLE PRECISION B(3,3), SIG_EQ(3,3), SIG_VIS(3,3)
      DOUBLE PRECISION FP(3,3), FE(3,3), INVFP(3,3)
      DOUBLE PRECISION BI1, JDET
      DOUBLE PRECISION STIFFN, P_PRESS, K_PRONY, H_PRESS
      INTEGER FLAG_VISC, HE_BASE, NET_BASE, N, NVISC, TAB
      INTEGER I, J
      DOUBLE PRECISION ZERO, ONE, SUM_STIFFN, EM20, THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, EM20=1.0D-20)
      PARAMETER (THIRD=ONE/3.0D0)
C
      DO I = 1, 3
        DO J = 1, 3
          SIG(I,J) = ZERO
        END DO
      END DO
C
      HE_BASE = 4
      IF (FLAG_HE .EQ. 1) THEN
        HE_BASE = HE_BASE + 14
      ELSE IF (FLAG_HE .EQ. 2) THEN
        HE_BASE = HE_BASE + 8
      ELSE IF (FLAG_HE .EQ. 3) THEN
        HE_BASE = HE_BASE + 6
      ELSE IF (FLAG_HE .EQ. 4) THEN
        HE_BASE = HE_BASE + 3
      ELSE IF (FLAG_HE .EQ. 5) THEN
        HE_BASE = HE_BASE + 9
      ELSE
        RETURN
      END IF
C
      SUM_STIFFN = ZERO
      TAB = HE_BASE
      DO N = 1, N_NETWORK
        STIFFN = ELCONLOC(TAB)
        FLAG_VISC = NINT(ELCONLOC(TAB+1))
        SUM_STIFFN = SUM_STIFFN + STIFFN
        IF (FLAG_VISC .EQ. 1) THEN
          NVISC = 5
        ELSE IF (FLAG_VISC .EQ. 2) THEN
          NVISC = 3
        ELSE IF (FLAG_VISC .EQ. 3) THEN
          NVISC = 3
        ELSE IF (FLAG_VISC .EQ. 4) THEN
          NVISC = 3
        ELSE
          TAB = TAB + 2
          CYCLE
        END IF
        TAB = TAB + 2 + NVISC
      END DO
C
      CALL PRODAAT(F, B)
      CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                  SIG_EQ, BI1, JDET)
      IF (N_NETWORK .GT. 0) THEN
        DO I = 1, 3
          DO J = 1, 3
            SIG(I,J) = (ONE - SUM_STIFFN) * SIG_EQ(I,J)
          END DO
        END DO
      ELSE
        DO I = 1, 3
          DO J = 1, 3
            SIG(I,J) = SIG_EQ(I,J)
          END DO
        END DO
      END IF
C
      TAB = HE_BASE
      DO N = 1, N_NETWORK
        STIFFN    = ELCONLOC(TAB)
        FLAG_VISC = NINT(ELCONLOC(TAB+1))
        IF (FLAG_VISC .EQ. 1) THEN
          NVISC = 5
        ELSE IF (FLAG_VISC .EQ. 2) THEN
          NVISC = 3
        ELSE IF (FLAG_VISC .EQ. 3) THEN
          NVISC = 3
        ELSE IF (FLAG_VISC .EQ. 4) THEN
          NVISC = 3
        ELSE
          TAB = TAB + 2
          CYCLE
        END IF
C
        NET_BASE = 2 + 12*(N-1)
        IF (FLAG_VISC .EQ. 4) THEN
          GOTO 900
        END IF
        FP(1,1) = STATEV(NET_BASE+1)
        FP(2,2) = STATEV(NET_BASE+2)
        FP(3,3) = STATEV(NET_BASE+3)
        FP(1,2) = STATEV(NET_BASE+4)
        FP(2,3) = STATEV(NET_BASE+5)
        FP(3,1) = STATEV(NET_BASE+6)
        FP(2,1) = STATEV(NET_BASE+7)
        FP(3,2) = STATEV(NET_BASE+8)
        FP(1,3) = STATEV(NET_BASE+9)
        CALL KMATINV3(FP, INVFP)
        CALL PRODMAT(F, INVFP, FE)
        CALL PRODAAT(FE, B)
        CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                    SIG_VIS, BI1, JDET)
        DO I = 1, 3
          DO J = 1, 3
            SIG(I,J) = SIG(I,J) + STIFFN * SIG_VIS(I,J)
          END DO
        END DO
        GOTO 901
C
  900   CONTINUE
        P_PRESS = THIRD*(SIG_EQ(1,1)+SIG_EQ(2,2)+SIG_EQ(3,3))
        K_PRONY = ELCONLOC(TAB+2)
        SIG_VIS(1,1) = (SIG_EQ(1,1)-P_PRESS) - STATEV(NET_BASE+1)
        SIG_VIS(2,2) = (SIG_EQ(2,2)-P_PRESS) - STATEV(NET_BASE+2)
        SIG_VIS(3,3) = (SIG_EQ(3,3)-P_PRESS) - STATEV(NET_BASE+3)
        SIG_VIS(1,2) = SIG_EQ(1,2) - STATEV(NET_BASE+4)
        SIG_VIS(2,3) = SIG_EQ(2,3) - STATEV(NET_BASE+5)
        SIG_VIS(3,1) = SIG_EQ(3,1) - STATEV(NET_BASE+6)
        H_PRESS = P_PRESS - STATEV(NET_BASE+7)
        SIG_VIS(2,1) = SIG_VIS(1,2)
        SIG_VIS(3,2) = SIG_VIS(2,3)
        SIG_VIS(1,3) = SIG_VIS(3,1)
        DO I = 1, 3
          DO J = 1, 3
            SIG(I,J) = SIG(I,J) + STIFFN * SIG_VIS(I,J)
          END DO
        END DO
        SIG(1,1) = SIG(1,1) + K_PRONY * H_PRESS
        SIG(2,2) = SIG(2,2) + K_PRONY * H_PRESS
        SIG(3,3) = SIG(3,3) + K_PRONY * H_PRESS
        SIG(1,1) = SIG(1,1) + STIFFN * P_PRESS
        SIG(2,2) = SIG(2,2) + STIFFN * P_PRESS
        SIG(3,3) = SIG(3,3) + STIFFN * P_PRESS
C
  901   CONTINUE
        TAB = TAB + 2 + NVISC
      END DO
C
      RETURN
      END


C ============================================================
C  PRF_UPDATE_OS — creep / Prony update (simplified for OptiStruct)
C  (calls same visc_models as CalculiX UMAT)
C ============================================================
      SUBROUTINE PRF_UPDATE_OS(F, DT, ELCONLOC, NTERMS,
     &                         N_NETWORK, FLAG_HE, FLAG_PL,
     &                         STATEV, NSTATV)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), DT, ELCONLOC(*), STATEV(*)
      INTEGER NTERMS, N_NETWORK, FLAG_HE, FLAG_PL, NSTATV
C
      DOUBLE PRECISION B(3,3), SIG(3,3), SIG_TRIAL(3,3)
      DOUBLE PRECISION FP(3,3), FP_OLD(3,3), FE(3,3)
      DOUBLE PRECISION INVFP(3,3), LB(3,3), FEDP(3,3)
      DOUBLE PRECISION DFP(3,3), SN(3,3), FP_NEW(3,3)
      DOUBLE PRECISION TRB, SB1, SB2, SB3, SB4, SB5, SB6
      DOUBLE PRECISION TBNORM, DGAMMA, FACTOR, AJ, BI1, JDET
      DOUBLE PRECISION STIFFN, A1, EXPC, EXPM, KSI, TAUREF
      DOUBLE PRECISION B0, EXPN, GAMMAOLD
      DOUBLE PRECISION G_I, TAU_I, ALPHA, DGAMMA_SUB, DGMAX, P_PRESS
      INTEGER FLAG_VISC, HE_BASE, NET_BASE, N, NVISC
      INTEGER TAB, I, J, NSUB, ISUB
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, EM20, THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (EM20=1.0D-20, THIRD=ONE/THREE)
      PARAMETER (DGMAX=1.0D-6)
C
      IF (N_NETWORK .LE. 0) RETURN
C
      HE_BASE = 4
      IF (FLAG_HE .EQ. 1) THEN
        HE_BASE = HE_BASE + 14
      ELSE IF (FLAG_HE .EQ. 2) THEN
        HE_BASE = HE_BASE + 8
      ELSE IF (FLAG_HE .EQ. 3) THEN
        HE_BASE = HE_BASE + 6
      ELSE IF (FLAG_HE .EQ. 4) THEN
        HE_BASE = HE_BASE + 3
      ELSE IF (FLAG_HE .EQ. 5) THEN
        HE_BASE = HE_BASE + 9
      ELSE
        RETURN
      END IF
C
      TAB = HE_BASE
      DO N = 1, N_NETWORK
        STIFFN    = ELCONLOC(TAB)
        FLAG_VISC = NINT(ELCONLOC(TAB+1))
C
        IF (FLAG_VISC .EQ. 1) THEN
          NVISC = 5
        ELSE IF (FLAG_VISC .EQ. 2) THEN
          NVISC = 3
        ELSE IF (FLAG_VISC .EQ. 3) THEN
          NVISC = 3
        ELSE IF (FLAG_VISC .EQ. 4) THEN
C         Prony: update h from current F
          CALL PRODAAT(F, B)
          CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                      SIG_TRIAL, BI1, JDET)
          TAU_I = ELCONLOC(TAB+3)
          ALPHA = EXP(-DT/MAX(EM20, TAU_I))
          NET_BASE = 2 + 12*(N-1)
          P_PRESS = THIRD * (SIG_TRIAL(1,1)+SIG_TRIAL(2,2)
     1                      +SIG_TRIAL(3,3))
          STATEV(NET_BASE+1) = ALPHA*STATEV(NET_BASE+1)
     1        + (ONE-ALPHA)*(SIG_TRIAL(1,1)-P_PRESS)
          STATEV(NET_BASE+2) = ALPHA*STATEV(NET_BASE+2)
     1        + (ONE-ALPHA)*(SIG_TRIAL(2,2)-P_PRESS)
          STATEV(NET_BASE+3) = ALPHA*STATEV(NET_BASE+3)
     1        + (ONE-ALPHA)*(SIG_TRIAL(3,3)-P_PRESS)
          STATEV(NET_BASE+4) = ALPHA*STATEV(NET_BASE+4)
     1        + (ONE-ALPHA)*SIG_TRIAL(1,2)
          STATEV(NET_BASE+5) = ALPHA*STATEV(NET_BASE+5)
     1        + (ONE-ALPHA)*SIG_TRIAL(2,3)
          STATEV(NET_BASE+6) = ALPHA*STATEV(NET_BASE+6)
     1        + (ONE-ALPHA)*SIG_TRIAL(3,1)
          STATEV(NET_BASE+7) = ALPHA*STATEV(NET_BASE+7)
     1        + (ONE-ALPHA)*P_PRESS
          TAB = TAB + 5
          CYCLE
        ELSE
          TAB = TAB + 2
          CYCLE
        END IF
C
C       BB / Sinh / Power: standard creep update
        NET_BASE = 2 + 12*(N-1)
        FP_OLD(1,1) = STATEV(NET_BASE+1)
        FP_OLD(2,2) = STATEV(NET_BASE+2)
        FP_OLD(3,3) = STATEV(NET_BASE+3)
        FP_OLD(1,2) = STATEV(NET_BASE+4)
        FP_OLD(2,3) = STATEV(NET_BASE+5)
        FP_OLD(3,1) = STATEV(NET_BASE+6)
        FP_OLD(2,1) = STATEV(NET_BASE+7)
        FP_OLD(3,2) = STATEV(NET_BASE+8)
        FP_OLD(1,3) = STATEV(NET_BASE+9)
C
        CALL KMATINV3(FP_OLD, INVFP)
        CALL PRODMAT(F, INVFP, FE)
        CALL PRODAAT(FE, B)
        CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                    SIG_TRIAL, BI1, JDET)
C
        SB1 = SIG_TRIAL(1,1)
        SB2 = SIG_TRIAL(2,2)
        SB3 = SIG_TRIAL(3,3)
        SB4 = SIG_TRIAL(1,2)
        SB5 = SIG_TRIAL(2,3)
        SB6 = SIG_TRIAL(3,1)
        TBNORM = SQRT(MAX(EM20,
     1     SB1*SB1 + SB2*SB2 + SB3*SB3
     2   + TWO*(SB4*SB4 + SB5*SB5 + SB6*SB6)))
C
        IF (FLAG_VISC .EQ. 1) THEN
          A1    = ELCONLOC(TAB+2)
          EXPC  = ELCONLOC(TAB+3)
          EXPM  = ELCONLOC(TAB+4)
          KSI   = ELCONLOC(TAB+5)
          TAUREF= ELCONLOC(TAB+6)
          CALL VISC_BB(FP_OLD, TBNORM, A1*DT, EXPC, EXPM, KSI,
     &                 TAUREF, DGAMMA)
        ELSE IF (FLAG_VISC .EQ. 2) THEN
          A1 = ELCONLOC(TAB+2)
          B0 = ELCONLOC(TAB+3)
          EXPN = ELCONLOC(TAB+4)
          CALL VISC_SINH(TBNORM, A1*DT, B0, EXPN, DGAMMA)
        ELSE
          A1    = ELCONLOC(TAB+2)
          EXPN  = ELCONLOC(TAB+3)
          EXPM  = ELCONLOC(TAB+4)
          GAMMAOLD = STATEV(NET_BASE+12)
          CALL VISC_POWER(TBNORM, A1*DT, EXPM, EXPN,
     &                    GAMMAOLD, DGAMMA)
        END IF
C
        NSUB = MAX(1, NINT(DGAMMA / DGMAX))
        DGAMMA_SUB = DGAMMA / DBLE(NSUB)
        DO ISUB = 1, NSUB
          FACTOR = DGAMMA_SUB / MAX(EM20, TBNORM)
          LB(1,1) = FACTOR * SB1
          LB(2,2) = FACTOR * SB2
          LB(3,3) = FACTOR * SB3
          LB(1,2) = FACTOR * SB4
          LB(2,3) = FACTOR * SB5
          LB(3,1) = FACTOR * SB6
          LB(2,1) = LB(1,2)
          LB(3,2) = LB(2,3)
          LB(1,3) = LB(3,1)
          CALL KMATINV3(FE, INVFP)
          CALL PRODMAT(LB, FE, FEDP)
          CALL PRODMAT(INVFP, FEDP, DFP)
          SN(1,1) = ONE + DFP(1,1)
          SN(2,2) = ONE + DFP(2,2)
          SN(3,3) = ONE + DFP(3,3)
          SN(1,2) = DFP(1,2)
          SN(2,3) = DFP(2,3)
          SN(3,1) = DFP(3,1)
          SN(2,1) = DFP(2,1)
          SN(3,2) = DFP(3,2)
          SN(1,3) = DFP(1,3)
          CALL PRODMAT(SN, FP_OLD, FP_NEW)
          FP_OLD(1,1) = FP_NEW(1,1)
          FP_OLD(2,2) = FP_NEW(2,2)
          FP_OLD(3,3) = FP_NEW(3,3)
          FP_OLD(1,2) = FP_NEW(1,2)
          FP_OLD(2,3) = FP_NEW(2,3)
          FP_OLD(3,1) = FP_NEW(3,1)
          FP_OLD(2,1) = FP_NEW(2,1)
          FP_OLD(3,2) = FP_NEW(3,2)
          FP_OLD(1,3) = FP_NEW(1,3)
          CALL KMATINV3(FP_OLD, INVFP)
          CALL PRODMAT(F, INVFP, FE)
        END DO
C
        STATEV(NET_BASE+1)  = FP_NEW(1,1)
        STATEV(NET_BASE+2)  = FP_NEW(2,2)
        STATEV(NET_BASE+3)  = FP_NEW(3,3)
        STATEV(NET_BASE+4)  = FP_NEW(1,2)
        STATEV(NET_BASE+5)  = FP_NEW(2,3)
        STATEV(NET_BASE+6)  = FP_NEW(3,1)
        STATEV(NET_BASE+7)  = FP_NEW(2,1)
        STATEV(NET_BASE+8)  = FP_NEW(3,2)
        STATEV(NET_BASE+9)  = FP_NEW(1,3)
        STATEV(NET_BASE+10) = DGAMMA
        STATEV(NET_BASE+11) = TBNORM
        IF (FLAG_VISC .EQ. 3) THEN
          STATEV(NET_BASE+12) = GAMMAOLD + DGAMMA
        END IF
        TAB = TAB + 2 + NVISC
      END DO
C
      RETURN
      END


C ============================================================
C  smatusr — Linear / small displacement (stub, not used for PRF)
C ============================================================
      SUBROUTINE smatusr(idu, nprop, prop, ndi, nshear, ntens,
     &                   smat, ieuid)
!DEC$ ATTRIBUTES DLLEXPORT :: smatusr
      IMPLICIT NONE
      INTEGER idu, nprop, ndi, nshear, ntens, ieuid
      DOUBLE PRECISION prop(nprop), smat(ntens,ntens)
C
C     PRF is nonlinear — smatusr should not be called.
C     Return zero stiffness matrix as fallback.
      INTEGER I, J
      DO I = 1, ntens
        DO J = 1, ntens
          SMAT(I,J) = 0.0D0
        END DO
      END DO
      RETURN
      END


C ============================================================
C  initusr — Initialize state variable labels (optional)
C ============================================================
      SUBROUTINE initusr(idu, nstate, cstate)
!DEC$ ATTRIBUTES DLLEXPORT :: initusr
      IMPLICIT NONE
      INTEGER idu, nstate
      CHARACTER*64 cstate(nstate)
      INTEGER I
C
      IF (nstate .GE. 1) cstate(1) = 'INIT_FLAG'
      IF (nstate .GE. 2) cstate(2) = 'J_DET'
      DO I = 3, MIN(nstate, 14)
        WRITE(cstate(I), '(A,I0)') 'NW1_STATEV_', I-2
      END DO
      DO I = 15, MIN(nstate, 26)
        WRITE(cstate(I), '(A,I0)') 'NW2_STATEV_', I-14
      END DO
      DO I = 27, MIN(nstate, 38)
        WRITE(cstate(I), '(A,I0)') 'NW3_STATEV_', I-26
      END DO
      RETURN
      END


C ============================================================
C  PRF_INIT_STATE_OS — initialize state variables
C  (Simplified version for OptiStruct)
C ============================================================
      SUBROUTINE PRF_INIT_STATE_OS(STATEV, N_NETWORK, FLAG_HE,
     &                              ELCONLOC)
      IMPLICIT NONE
      DOUBLE PRECISION STATEV(*), ELCONLOC(*)
      INTEGER N_NETWORK, FLAG_HE
      INTEGER N, BASE, TAB, HE_BASE, FLAG_VISC
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      STATEV(1) = -1.0D0
      STATEV(2) = ONE
C
      HE_BASE = 4
      IF (FLAG_HE .EQ. 1) THEN
        HE_BASE = HE_BASE + 14
      ELSE IF (FLAG_HE .EQ. 2) THEN
        HE_BASE = HE_BASE + 8
      ELSE IF (FLAG_HE .EQ. 3) THEN
        HE_BASE = HE_BASE + 6
      ELSE IF (FLAG_HE .EQ. 4) THEN
        HE_BASE = HE_BASE + 3
      ELSE IF (FLAG_HE .EQ. 5) THEN
        HE_BASE = HE_BASE + 9
      END IF
C
      TAB = HE_BASE
      DO N = 1, N_NETWORK
        BASE = 2 + 12*(N-1)
        FLAG_VISC = NINT(ELCONLOC(TAB+1))
        IF (FLAG_VISC .EQ. 4) THEN
          STATEV(BASE+1) = ZERO
          STATEV(BASE+2) = ZERO
          STATEV(BASE+3) = ZERO
          STATEV(BASE+4) = ZERO
          STATEV(BASE+5) = ZERO
          STATEV(BASE+6) = ZERO
          STATEV(BASE+7) = ZERO
          TAB = TAB + 5
        ELSE
          STATEV(BASE+1) = ONE
          STATEV(BASE+2) = ONE
          STATEV(BASE+3) = ONE
          STATEV(BASE+4) = ZERO
          STATEV(BASE+5) = ZERO
          STATEV(BASE+6) = ZERO
          STATEV(BASE+7) = ZERO
          STATEV(BASE+8) = ZERO
          STATEV(BASE+9) = ZERO
          STATEV(BASE+10)= ZERO
          STATEV(BASE+11)= ZERO
          STATEV(BASE+12)= 1.0D-20
          TAB = TAB + 5
        END IF
      END DO
      RETURN
      END
