C ============================================================
C  WHTOOLs MATCALIB 2026 — Parallel Rheological Framework UMAT
C  for CalculiX ccx_custom
C
C  Implements the Multi-Network Foam (MNF) / PRF material
C  model from OpenRadioss LAW100 (sigeps100.F90).
C
C  Ported sub-models:
C    Polynomial hyperelastic (sigpoly.F)    via polynomial_stress.f
C    Arruda-Boyce 8-chain (sigaboyce.F)     via arruda_boyce.f
C    Bergstrom-Boyce creep (viscbb.F)       via visc_models.f
C    Power law creep (viscpower.F)          via visc_models.f
C    Sinh creep (viscsinh.F)                via visc_models.f
C    Prony series (linear viscoelastic)     via visc_models.f
C
C  Keyword input (CalculiX *USER MATERIAL):
C    *USER MATERIAL, CONSTANTS=NC
C    N_NETWORK, FLAG_HE, FLAG_PL
C    [hyperelastic constants (see below)]
C    [viscous network constants (repeated N_NETWORK times)]
C    G, RBULK
C
C  where:
C    N_NETWORK = number of viscous networks (integer, ≥0)
C    FLAG_HE   = 1 (Polynomial), 2 (Arruda-Boyce)
C    FLAG_PL   = 0 (no plasticity) — future use
C
C  Hyperelastic constants (after N_NETWORK, FLAG_HE, FLAG_PL):
C
C    FLAG_HE=1 (Polynomial, 14 constants):
C      C10, C01, C20, C11, C02, C30, C21, C12, C03, D1, D2, D3,
C      IHYPER, RBULK_BACKUP
C
C    FLAG_HE=2 (Arruda-Boyce, 8 constants):
C      C1, C2, C3, C4, C5, MU, D, BETA
C
C  Per viscous network:
C    STIFFN, FLAG_VISC, [visc parameters...]
C
C    FLAG_VISC=1 (Bergstrom-Boyce, 5 params):
C      A1, EXPC, EXPM, KSI, TAUREF
C    FLAG_VISC=2 (Sinh, 3 params):
C      A1, B0, EXPN
C    FLAG_VISC=3 (Power law, 3 params):
C      A1, EXPN, EXPM
C
C  Final:
C    G (shear modulus for wave speed), RBULK (bulk modulus)
C
C  State variables (*DEPVAR):
C    NV = 2 + 12 * N_NETWORK
C    statev(1) = initialised flag
C    statev(2) = J (determinant, stored)
C    Per viscous network (FLAG_VISC=1/2/3):
C      statev(2+12*n+1..9)  = Fv(3,3) viscous def grad
C      statev(2+12*n+10)    = dgamma (creep increment)
C      statev(2+12*n+11)    = tbnorm (effective stress norm)
C      statev(2+12*n+12)    = gamma_old (cumulative, Power law)
C    Per Prony term (FLAG_VISC=4):
C      statev(2+12*n+1..6)  = h(1..6) hereditary stress tensor
C ============================================================
      SUBROUTINE umat_user(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xokl,voj,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mi,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,pnewdt,ipkon)
C
      IMPLICIT NONE
C
      CHARACTER*80 amat
      INTEGER ithermal(*),icmd,kode,ielas,iel,iint,nstate_,mi(*),
     &        iorien,ipkon(*)
      REAL*8 elconloc(*),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &       vj,t1l,dtime,xkl(3,3),xokl(3,3),voj,pgauss(3),orab(7,*),
     &       time,ttime,pnewdt
      REAL*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
C
C ---- local variables -----------------------------------------
      DOUBLE PRECISION F(3,3), SIG_CAUCHY(3,3), PK2(3,3)
      DOUBLE PRECISION FINV(3,3), FINV_T(3,3), VEC(3,3)
      DOUBLE PRECISION AJ, DET
      DOUBLE PRECISION S_REF(6), S_PERT(6)
      DOUBLE PRECISION FPERT(3,3), DELTAF(3,3)
      DOUBLE PRECISION STATEV(100)     ! working copy (max 100)
      DOUBLE PRECISION EPS, IWK
      INTEGER NTERMS, NTERMS_BASE
      INTEGER N_NETWORK, FLAG_HE, FLAG_PL
      INTEGER NHE, I, J, ICOL, K, P, Q, NSTATV
      INTEGER NPOLY, ISTART, IP
      PARAMETER (EPS=1.0D-8)
C
C     maximum state variables per point (must match *DEPVAR)
      INTEGER MAXSV
      PARAMETER (MAXSV=100)
C
C     Voigt index mapping: 1→11, 2→22, 3→33, 4→12, 5→13, 6→23
      INTEGER VPMAP(2,6)
      DATA VPMAP /1,1, 2,2, 3,3, 1,2, 1,3, 2,3/
C
C ---- 1. determine number of constants ---------------------------
      NTERMS = -kode - 100
      IF (NTERMS .LT. 5) THEN
        WRITE(*,*) 'ERROR umat_user: need at least 5 constants, got ',
     &             NTERMS
        pnewdt = 0.0D0
        RETURN
      END IF
C
C ---- 2. parse top-level parameters ------------------------------
      N_NETWORK = NINT(elconloc(1))
      FLAG_HE   = NINT(elconloc(2))
      FLAG_PL   = NINT(elconloc(3))
C
C     compute number of state variables
      NSTATV = 2 + 12 * N_NETWORK
      IF (NSTATV .GT. MAXSV) THEN
        WRITE(*,*) 'ERROR umat_user: NSTATV=', NSTATV,
     &             ' exceeds MAXSV=', MAXSV
        pnewdt = 0.0D0
        RETURN
      END IF
C
C ---- 3. copy deformation gradient (CalculiX stores F transposed) -
      DO I = 1, 3
        DO J = 1, 3
          F(I,J) = xkl(J,I)
        END DO
      END DO
C
C ---- 4. read / initialise state variables -----------------------
      DO I = 1, NSTATV
        STATEV(I) = xstate(I, 1, 1)
      END DO
C
      IF (time .EQ. 0.0D0 .OR. STATEV(1) .EQ. 0.0D0) THEN
        CALL PRF_INIT_STATE(STATEV, N_NETWORK, NSTATV, FLAG_HE,
     &                       ELCONLOC)
      END IF
C
C ---- 5. compute J = det(F) --------------------------------------
      AJ = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
     1   - F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
     2   + F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
      STATEV(2) = AJ
C
C ---- 6. update internal variables (creep update) ----------------
      CALL PRF_UPDATE(F, DTIME, ELCONLOC, NTERMS, N_NETWORK,
     &                FLAG_HE, FLAG_PL, STATEV, NSTATV)
C
C ---- 7. compute total Cauchy stress (with updated state) -------
      CALL PRF_EVAL(F, DTIME, ELCONLOC, NTERMS, N_NETWORK, FLAG_HE,
     &              FLAG_PL, STATEV, NSTATV, SIG_CAUCHY)
C
C ---- 8. Cauchy → PK2: S = J * F^-1 * sigma * F^-T ------------
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
      stre(1) = PK2(1,1)
      stre(2) = PK2(2,2)
      stre(3) = PK2(3,3)
      stre(4) = PK2(1,2)
      stre(5) = PK2(1,3)
      stre(6) = PK2(2,3)
C
C ---- 9. save reference PK2 for tangent ------------------------
      S_REF(1) = PK2(1,1)
      S_REF(2) = PK2(2,2)
      S_REF(3) = PK2(3,3)
      S_REF(4) = PK2(1,2)
      S_REF(5) = PK2(1,3)
      S_REF(6) = PK2(2,3)
C
C ---- 11. numerical tangent stiffness ----------------------------
      IF(icmd .NE. 3) THEN
        DO I = 1, 21
          stiff(I) = 0.0D0
        END DO
C
C       recompute FINV_T for perturbation
        DO I = 1, 3
          DO J = 1, 3
            FINV_T(I,J) = FINV(J,I)
          END DO
        END DO
C
        DO ICOL = 1, 6
          P = VPMAP(1, ICOL)
          Q = VPMAP(2, ICOL)
C
C         perturbation direction: δF = eps * F^-T * (e_p ⊗ e_q)
          DO I = 1, 3
            DO J = 1, 3
              DELTAF(I,J) = 0.0D0
            END DO
          END DO
C         Miehe numerical tangent perturbation:
C         δF = eps * F^{-T} * (e_p ⊗ e_q)  →  δF_{i,q} = eps * (F^{-T})_{i,p}
          IF(P .EQ. Q) THEN
            DO I = 1, 3
              DELTAF(I,P) = EPS * FINV_T(I,P)
            END DO
          ELSE
            DO I = 1, 3
              DELTAF(I,Q) = EPS * FINV_T(I,P)
            END DO
          END IF
C
C         perturbed F
          DO I = 1, 3
            DO J = 1, 3
              FPERT(I,J) = F(I,J) + DELTAF(I,J)
            END DO
          END DO
C
C         frozen tangent: PRF_EVAL with converged Fv
          CALL PRF_EVAL(FPERT, DTIME, ELCONLOC, NTERMS, N_NETWORK,
     &                  FLAG_HE, FLAG_PL, STATEV, NSTATV,
     &                  SIG_CAUCHY)
C
C         Cauchy → PK2 for perturbed state
          CALL INV3X3(FPERT, FINV, DET)
          DO I = 1, 3
            DO J = 1, 3
              FINV_T(I,J) = FINV(J,I)
            END DO
          END DO
          DET = FPERT(1,1)*(FPERT(2,2)*FPERT(3,3)-FPERT(2,3)*FPERT(3,2))
     1        - FPERT(1,2)*(FPERT(2,1)*FPERT(3,3)-FPERT(2,3)*FPERT(3,1))
     2        + FPERT(1,3)*(FPERT(2,1)*FPERT(3,2)-FPERT(2,2)*FPERT(3,1))
          CALL MATMUL_3X3(FINV, SIG_CAUCHY, VEC)
          CALL MATMUL_3X3(VEC, FINV_T, PK2)
          DO I = 1, 3
            DO J = 1, 3
              PK2(I,J) = DET * PK2(I,J)
            END DO
          END DO
C
          S_PERT(1) = PK2(1,1)
          S_PERT(2) = PK2(2,2)
          S_PERT(3) = PK2(3,3)
          S_PERT(4) = PK2(1,2)
          S_PERT(5) = PK2(1,3)
          S_PERT(6) = PK2(2,3)
C
C         finite difference: dS_ij / dE_pq
          DO I = 1, 6
            IWK = (S_PERT(I) - S_REF(I)) / EPS
            IF(ICOL .EQ. 1) THEN
              IF(I .EQ. 1) stiff(1)  = IWK
              IF(I .EQ. 2) stiff(2)  = IWK
              IF(I .EQ. 3) stiff(4)  = IWK
              IF(I .EQ. 4) stiff(7)  = IWK
              IF(I .EQ. 5) stiff(11) = IWK
              IF(I .EQ. 6) stiff(16) = IWK
            ELSE IF(ICOL .EQ. 2) THEN
              IF(I .EQ. 1) stiff(2)  = stiff(2)  + IWK
              IF(I .EQ. 2) stiff(3)  = IWK
              IF(I .EQ. 3) stiff(5)  = IWK
              IF(I .EQ. 4) stiff(8)  = IWK
              IF(I .EQ. 5) stiff(12) = IWK
              IF(I .EQ. 6) stiff(17) = IWK
            ELSE IF(ICOL .EQ. 3) THEN
              IF(I .EQ. 1) stiff(4)  = stiff(4)  + IWK
              IF(I .EQ. 2) stiff(5)  = stiff(5)  + IWK
              IF(I .EQ. 3) stiff(6)  = IWK
              IF(I .EQ. 4) stiff(9)  = IWK
              IF(I .EQ. 5) stiff(13) = IWK
              IF(I .EQ. 6) stiff(18) = IWK
            ELSE IF(ICOL .EQ. 4) THEN
              IF(I .EQ. 1) stiff(7)  = stiff(7)  + IWK
              IF(I .EQ. 2) stiff(8)  = stiff(8)  + IWK
              IF(I .EQ. 3) stiff(9)  = stiff(9)  + IWK
              IF(I .EQ. 4) stiff(10) = IWK
              IF(I .EQ. 5) stiff(14) = IWK
              IF(I .EQ. 6) stiff(19) = IWK
            ELSE IF(ICOL .EQ. 5) THEN
              IF(I .EQ. 1) stiff(11) = stiff(11) + IWK
              IF(I .EQ. 2) stiff(12) = stiff(12) + IWK
              IF(I .EQ. 3) stiff(13) = stiff(13) + IWK
              IF(I .EQ. 4) stiff(14) = stiff(14) + IWK
              IF(I .EQ. 5) stiff(15) = IWK
              IF(I .EQ. 6) stiff(20) = IWK
            ELSE IF(ICOL .EQ. 6) THEN
              IF(I .EQ. 1) stiff(16) = stiff(16) + IWK
              IF(I .EQ. 2) stiff(17) = stiff(17) + IWK
              IF(I .EQ. 3) stiff(18) = stiff(18) + IWK
              IF(I .EQ. 4) stiff(19) = stiff(19) + IWK
              IF(I .EQ. 5) stiff(20) = stiff(20) + IWK
              IF(I .EQ. 6) stiff(21) = IWK
            END IF
          END DO
        END DO
C
C       symmetrize off-diagonal terms
        stiff(2)  = stiff(2)  / 2.0D0
        stiff(4)  = stiff(4)  / 2.0D0
        stiff(5)  = stiff(5)  / 2.0D0
        stiff(7)  = stiff(7)  / 2.0D0
        stiff(8)  = stiff(8)  / 2.0D0
        stiff(9)  = stiff(9)  / 2.0D0
        stiff(11) = stiff(11) / 2.0D0
        stiff(12) = stiff(12) / 2.0D0
        stiff(13) = stiff(13) / 2.0D0
        stiff(14) = stiff(14) / 2.0D0
        stiff(16) = stiff(16) / 2.0D0
        stiff(17) = stiff(17) / 2.0D0
        stiff(18) = stiff(18) / 2.0D0
        stiff(19) = stiff(19) / 2.0D0
        stiff(20) = stiff(20) / 2.0D0
      END IF
C
C ---- 12. store state variables back ----------------------------
      DO I = 1, NSTATV
        xstate(I, 1, 1) = STATEV(I)
      END DO
C
      pnewdt = -1.0D0
      RETURN
      END


C ============================================================
C  PRF_INIT_STATE - initialize state variables
C ============================================================
      SUBROUTINE PRF_INIT_STATE(STATEV, N_NETWORK, NSTATV, FLAG_HE,
     &                           ELCONLOC)
      IMPLICIT NONE
      DOUBLE PRECISION STATEV(*), ELCONLOC(*)
      INTEGER N_NETWORK, NSTATV, FLAG_HE
      INTEGER N, BASE, TAB, HE_BASE, FLAG_VISC
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C     statev(1) = init flag (set to -1 to indicate initialised)
C     statev(2) = J (will be set in main UMAT)
      STATEV(1) = -1.0D0
      STATEV(2) = ONE
C
C     find HE_BASE (same as PRF_EVAL)
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
C     initialise state variables per network
      TAB = HE_BASE
      DO N = 1, N_NETWORK
        BASE = 2 + 12*(N-1)
        FLAG_VISC = NINT(ELCONLOC(TAB+1))
C       Prony (FLAG_VISC=4): h_i starts at 0
C       BB/Sinh/Power (FLAG_VISC=1/2/3): Fv starts at I
        IF (FLAG_VISC .EQ. 4) THEN
          STATEV(BASE+1) = ZERO  ! h_11
          STATEV(BASE+2) = ZERO  ! h_22
          STATEV(BASE+3) = ZERO  ! h_33
          STATEV(BASE+4) = ZERO  ! h_12
          STATEV(BASE+5) = ZERO  ! h_23
          STATEV(BASE+6) = ZERO  ! h_31
          TAB = TAB + 5           ! 2 header + 3 params
        ELSE
          STATEV(BASE+1) = ONE   ! Fv_11
          STATEV(BASE+2) = ONE   ! Fv_22
          STATEV(BASE+3) = ONE   ! Fv_33
          STATEV(BASE+4) = ZERO  ! Fv_12
          STATEV(BASE+5) = ZERO  ! Fv_23
          STATEV(BASE+6) = ZERO  ! Fv_31
          STATEV(BASE+7) = ZERO  ! Fv_21
          STATEV(BASE+8) = ZERO  ! Fv_32
          STATEV(BASE+9) = ZERO  ! Fv_13
          STATEV(BASE+10)= ZERO  ! dgamma
          STATEV(BASE+11)= ZERO  ! tbnorm
          STATEV(BASE+12)= 1.0D-20 ! gamma_old
          TAB = TAB + 2 + 3      ! default advance (will be overwritten)
        END IF
      END DO
C
      RETURN
      END


C ============================================================
C  PRF_UPDATE - creep update for all viscous networks
C  Modifies STATEV in-place (updates Fv, dgamma, tbnorm)
C ============================================================
      SUBROUTINE PRF_UPDATE(F, DT, ELCONLOC, NTERMS,
     &                      N_NETWORK, FLAG_HE, FLAG_PL,
     &                      STATEV, NSTATV)
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
       DOUBLE PRECISION G_I, TAU_I, ALPHA, DGAMMA_SUB, P_PRESS
       INTEGER FLAG_VISC, HE_BASE, NET_BASE, N, NVISC
       INTEGER TAB, I, J, NSUB, ISUB
        DOUBLE PRECISION ZERO, ONE, TWO, THREE, EM20, THIRD, DGMAX
        PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
       PARAMETER (EM20=1.0D-20, THIRD=ONE/THREE)
       PARAMETER (DGMAX=1.0D-6)
C
C     skip if no viscous networks
      IF (N_NETWORK .LE. 0) RETURN
C
C     find hyperelastic parameter start
      HE_BASE = 4  ! after N_NETWORK, FLAG_HE, FLAG_PL
      IF (FLAG_HE .EQ. 1) THEN
        HE_BASE = HE_BASE + 14  ! polynomial
      ELSE IF (FLAG_HE .EQ. 2) THEN
        HE_BASE = HE_BASE + 8   ! Arruda-Boyce
      ELSE IF (FLAG_HE .EQ. 3) THEN
        HE_BASE = HE_BASE + 6   ! Yeoh
      ELSE IF (FLAG_HE .EQ. 4) THEN
        HE_BASE = HE_BASE + 3   ! Gent
      ELSE IF (FLAG_HE .EQ. 5) THEN
        HE_BASE = HE_BASE + 9   ! Ogden
      ELSE
        RETURN  ! unknown flag_he
      END IF
C
      IF (FLAG_PL .EQ. 1) THEN
        HE_BASE = HE_BASE + 5   ! plasticity params
      END IF
C
C     loop over viscous networks
      TAB = HE_BASE
      DO N = 1, N_NETWORK
        STIFFN    = ELCONLOC(TAB)
        FLAG_VISC = NINT(ELCONLOC(TAB+1))
C
C       determine NVISC (number of viscous params + 2 header)
        IF (FLAG_VISC .EQ. 1) THEN
          NVISC = 5   ! BB: A1, EXPC, EXPM, KSI, TAUREF
        ELSE IF (FLAG_VISC .EQ. 2) THEN
          NVISC = 3   ! Sinh: A1, B0, EXPN
        ELSE IF (FLAG_VISC .EQ. 3) THEN
          NVISC = 3   ! Power: A1, EXPN, EXPM
        ELSE IF (FLAG_VISC .EQ. 4) THEN
C         Prony: update h from trial stress
          CALL PRODAAT(F, B)
          CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                      SIG_TRIAL, BI1, JDET)
          TAU_I = ELCONLOC(TAB+3)
          ALPHA = EXP(-DT/MAX(EM20, TAU_I))
          P_PRESS = THIRD * (SIG_TRIAL(1,1)+SIG_TRIAL(2,2)+SIG_TRIAL(3,3))
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
        ELSE
          TAB = TAB + 2
          CYCLE
        END IF
C
C       extract creep parameters
        A1 = ELCONLOC(TAB+2)
        IF (FLAG_VISC .EQ. 1) THEN
          EXPC   = ELCONLOC(TAB+3)
          EXPM   = ELCONLOC(TAB+4)
          KSI    = ELCONLOC(TAB+5)
          TAUREF = ELCONLOC(TAB+6)
        ELSE IF (FLAG_VISC .EQ. 2) THEN
          B0  = ELCONLOC(TAB+3)
          EXPN = ELCONLOC(TAB+4)
        ELSE IF (FLAG_VISC .EQ. 3) THEN
          EXPN = ELCONLOC(TAB+3)
          EXPM = ELCONLOC(TAB+4)
        END IF
C
C       read Fv from state variables
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
C       trial elastic b: fe = F * inv(Fp), b = fe * fe^T
        CALL KMATINV3(FP_OLD, INVFP)
        CALL PRODMAT(F, INVFP, FE)
        CALL PRODAAT(FE, B)
C
C       compute trial Cauchy stress (hyperelastic)
        CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                    SIG_TRIAL, BI1, JDET)
C
C       scale by stiffn
        DO I = 1, 3
          DO J = 1, 3
            SIG_TRIAL(I,J) = STIFFN * SIG_TRIAL(I,J)
          END DO
        END DO
C
C       deviatoric stress norm
        TRB = THIRD * (SIG_TRIAL(1,1)+SIG_TRIAL(2,2)+SIG_TRIAL(3,3))
        SB1 = SIG_TRIAL(1,1) - TRB
        SB2 = SIG_TRIAL(2,2) - TRB
        SB3 = SIG_TRIAL(3,3) - TRB
        SB4 = SIG_TRIAL(1,2)
        SB5 = SIG_TRIAL(2,3)
        SB6 = SIG_TRIAL(3,1)
        TBNORM = SQRT(MAX(EM20,
     1     SB1*SB1 + SB2*SB2 + SB3*SB3
     2   + TWO*(SB4*SB4 + SB5*SB5 + SB6*SB6)))
C
C       compute creep increment
        IF (FLAG_VISC .EQ. 1) THEN
C         Bergstrom-Boyce
          CALL VISC_BB(FP_OLD, TBNORM, A1*DT, EXPC, EXPM, KSI,
     &                 TAUREF, DGAMMA)
        ELSE IF (FLAG_VISC .EQ. 2) THEN
C         Sinh
          CALL VISC_SINH(TBNORM, A1*DT, B0, EXPN, DGAMMA)
        ELSE IF (FLAG_VISC .EQ. 3) THEN
C         Power law
          GAMMAOLD = STATEV(NET_BASE+12)
          CALL VISC_POWER(TBNORM, A1*DT, EXPM, EXPN,
     &                    GAMMAOLD, DGAMMA)
        END IF
C
C       sub-stepping for explicit Fv update (prevents overshoot)
        NSUB = MAX(1, NINT(DGAMMA / DGMAX))
        DGAMMA_SUB = DGAMMA / DBLE(NSUB)
        DO ISUB = 1, NSUB
C
C         flow direction: Lv_sub = (dg_sub / tau_norm) * dev(sigma)
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
C
C         dfp = inv(Fe) * Lv_sub * Fe
          CALL KMATINV3(FE, INVFP)
          CALL PRODMAT(LB, FE, FEDP)
          CALL PRODMAT(INVFP, FEDP, DFP)
C
C         Fv_new = (I + dfp) * Fv_old
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
C
C         update Fv_old and Fe for next sub-step
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
C       store updated state
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
C
        TAB = TAB + 2 + NVISC
      END DO
C
      RETURN
      END


C ============================================================
C  PRF_EVAL - compute total Cauchy stress (frozen internal vars)
C  Does NOT modify STATEV.
C ============================================================
      SUBROUTINE PRF_EVAL(F, DT, ELCONLOC, NTERMS,
     &                    N_NETWORK, FLAG_HE, FLAG_PL,
     &                    STATEV, NSTATV, SIG)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), DT, ELCONLOC(*), STATEV(*), SIG(3,3)
      INTEGER NTERMS, N_NETWORK, FLAG_HE, FLAG_PL, NSTATV
C
      DOUBLE PRECISION B(3,3), SIG_EQ(3,3), SIG_VIS(3,3)
      DOUBLE PRECISION FP(3,3), FE(3,3), INVFP(3,3)
      DOUBLE PRECISION BI1, JDET
      DOUBLE PRECISION STIFFN
      DOUBLE PRECISION P_PRESS, K_PRONY, H_PRESS, ALPHA, TAU_I
      INTEGER FLAG_VISC, HE_BASE, NET_BASE, N, NVISC, TAB
      INTEGER I, J
      DOUBLE PRECISION ZERO, ONE, SUM_STIFFN, EM20, THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, EM20=1.0D-20)
      PARAMETER (THIRD=1.0D0/3.0D0)
C
C     zero total stress
      DO I = 1, 3
        DO J = 1, 3
          SIG(I,J) = ZERO
        END DO
      END DO
C
C     find hyperelastic parameter start
      HE_BASE = 4
      IF (FLAG_HE .EQ. 1) THEN
        HE_BASE = 4 + 14
      ELSE IF (FLAG_HE .EQ. 2) THEN
        HE_BASE = 4 + 8
      ELSE IF (FLAG_HE .EQ. 3) THEN
        HE_BASE = 4 + 6
      ELSE IF (FLAG_HE .EQ. 4) THEN
        HE_BASE = 4 + 3
      ELSE IF (FLAG_HE .EQ. 5) THEN
        HE_BASE = 4 + 9
      ELSE
        RETURN
      END IF
      IF (FLAG_PL .EQ. 1) THEN
        HE_BASE = HE_BASE + 5
      END IF
C
C     ---- sum stiffn for equilibrium scaling ----
C     PRF: total = (1 - sum(stiffn)) * sigma_eq + sum(stiffn_i * sigma_vis_i)
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
          NVISC = 3   ! Prony: g, k, tau
        ELSE
          TAB = TAB + 2
          CYCLE
        END IF
        TAB = TAB + 2 + NVISC
      END DO
C
C     ---- equilibrium network (full b = F*F^T), scaled by (1-sum(stiffn)) ----
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
C       no viscous networks → pure hyperelastic
        DO I = 1, 3
          DO J = 1, 3
            SIG(I,J) = SIG_EQ(I,J)
          END DO
        END DO
      END IF
C
C     ---- viscous networks (b_e = Fe * Fe^T, Fe = F * inv(Fv)) ----
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
          NVISC = 3   ! Prony
        ELSE
          TAB = TAB + 2
          CYCLE
        END IF
C
C       read state variables
        NET_BASE = 2 + 12*(N-1)
        IF (FLAG_VISC .EQ. 4) THEN
C         Prony: skip Fv/Fe computation, use stored h directly
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
C
C       compute elastic left Cauchy-Green be = Fe * Fe^T
        CALL KMATINV3(FP, INVFP)
        CALL PRODMAT(F, INVFP, FE)
        CALL PRODAAT(FE, B)
C
C       hyperelastic stress contribution
  900   CONTINUE
        IF (FLAG_VISC .EQ. 4) THEN
C       Prony: evaluate stress using stored h, then update h with current SIG_EQ
C       After init: SIG = (1-Sg)*SIG_EQ + Sg*(dev(SIG_EQ)-h_dev+vol) = SIG_EQ - Sg*h_dev
        P_PRESS = THIRD*(SIG_EQ(1,1)+SIG_EQ(2,2)+SIG_EQ(3,3))
        K_PRONY = ELCONLOC(TAB+2)
        TAU_I   = ELCONLOC(TAB+3)
        ALPHA   = EXP(-DT/MAX(EM20, TAU_I))
C       Compute viscoelastic contribution from stored h
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
C       Compensate for (1-g)*SIG_EQ initial scaling of volumetric part:
C       Prony g_i only relaxes deviatoric, not volumetric. Add back g*p.
        SIG(1,1) = SIG(1,1) + STIFFN * P_PRESS
        SIG(2,2) = SIG(2,2) + STIFFN * P_PRESS
        SIG(3,3) = SIG(3,3) + STIFFN * P_PRESS
        ELSE
          CALL HYPER_STRESS(FLAG_HE, B, ELCONLOC(4), NTERMS,
     &                      SIG_VIS, BI1, JDET)
          DO I = 1, 3
            DO J = 1, 3
              SIG(I,J) = SIG(I,J) + STIFFN * SIG_VIS(I,J)
            END DO
          END DO
        END IF
        TAB = TAB + 2 + NVISC
      END DO
C
      RETURN
      END


C ============================================================
C  HYPER_STRESS - hyperelastic stress dispatcher
C  FLAG_HE=1 → POLY_STRESS (polynomial)
C  FLAG_HE=2 → ARRUDA_BOYCE (Arruda-Boyce 8-chain)
C  FLAG_HE=3 → YEOH_STRESS (Yeoh / Reduced Polynomial)
C  FLAG_HE=4 → GENT_STRESS (Gent locking model)
C  FLAG_HE=5 → OGDEN_STRESS (Ogden N=3)
C ============================================================
      SUBROUTINE HYPER_STRESS(FLAG_HE, B, ELC_BASE, NTERMS,
     &                        SIG, BI1, JDET)
      IMPLICIT NONE
      INTEGER FLAG_HE, NTERMS
      DOUBLE PRECISION B(3,3), ELC_BASE(*), SIG(3,3), BI1, JDET
C
      DOUBLE PRECISION C10, C01, C20, C11, C02, C30, C21, C12, C03
      DOUBLE PRECISION D1, D2, D3, RBULK, RIHYPER
      DOUBLE PRECISION C1, C2, C3, C4, C5, MU, D, BETA
      DOUBLE PRECISION BI2
      INTEGER IHYPER
C
      IF (FLAG_HE .EQ. 1) THEN
C       Polynomial hyperelastic
        C10 = ELC_BASE(1)
        C01 = ELC_BASE(2)
        C20 = ELC_BASE(3)
        C11 = ELC_BASE(4)
        C02 = ELC_BASE(5)
        C30 = ELC_BASE(6)
        C21 = ELC_BASE(7)
        C12 = ELC_BASE(8)
        C03 = ELC_BASE(9)
        D1  = ELC_BASE(10)
        D2  = ELC_BASE(11)
        D3  = ELC_BASE(12)
        RIHYPER = ELC_BASE(13)
        IHYPER = NINT(RIHYPER)
        RBULK = ELC_BASE(14)
        CALL POLY_STRESS(B, C10, C01, C20, C11, C02, C30, C21, C12,
     &                   C03, D1, D2, D3, SIG, BI1, BI2, JDET,
     &                   IHYPER, RBULK)
      ELSE IF (FLAG_HE .EQ. 2) THEN
C       Arruda-Boyce 8-chain
        C1   = ELC_BASE(1)
        C2   = ELC_BASE(2)
        C3   = ELC_BASE(3)
        C4   = ELC_BASE(4)
        C5   = ELC_BASE(5)
        MU   = ELC_BASE(6)
        D    = ELC_BASE(7)
        BETA = ELC_BASE(8)
        CALL ARRUDA_BOYCE(B, C1, C2, C3, C4, C5, MU, D, BETA,
     &                    SIG, BI1, JDET)
      ELSE IF (FLAG_HE .EQ. 3) THEN
C       Yeoh (Reduced Polynomial 3rd order)
        CALL YEOH_STRESS(B,
     1       ELC_BASE(1), ELC_BASE(2), ELC_BASE(3),
     2       ELC_BASE(4), ELC_BASE(5), ELC_BASE(6),
     3       SIG, BI1, JDET)
      ELSE IF (FLAG_HE .EQ. 4) THEN
C       Gent (locking model)
        CALL GENT_STRESS(B,
     1       ELC_BASE(1), ELC_BASE(2), ELC_BASE(3),
     2       SIG, BI1, JDET)
      ELSE IF (FLAG_HE .EQ. 5) THEN
C       Ogden (N=3 terms)
        CALL OGDEN_STRESS(B,
     1       ELC_BASE(1), ELC_BASE(2),
     2       ELC_BASE(3), ELC_BASE(4),
     3       ELC_BASE(5), ELC_BASE(6),
     4       ELC_BASE(7), ELC_BASE(8), ELC_BASE(9),
     5       SIG, BI1, JDET)
      END IF
C
      RETURN
      END


C ============================================================
C  3×3 matrix utilities
C ============================================================
      SUBROUTINE INV3X3(A, AINV, DET)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), AINV(3,3), DET
C
      DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1    - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2    + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      IF(ABS(DET) .LT. 1.0D-20) DET = 1.0D-20
C
      AINV(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
      AINV(2,1) = (A(3,2)*A(1,3)-A(1,2)*A(3,3))/DET
      AINV(3,1) = (A(1,2)*A(2,3)-A(2,2)*A(1,3))/DET
      AINV(1,2) = (A(2,3)*A(3,1)-A(2,1)*A(3,3))/DET
      AINV(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
      AINV(3,2) = (A(2,1)*A(1,3)-A(1,1)*A(2,3))/DET
      AINV(1,3) = (A(2,1)*A(3,2)-A(3,1)*A(2,2))/DET
      AINV(2,3) = (A(3,1)*A(1,2)-A(1,1)*A(3,2))/DET
      AINV(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
C
      RETURN
      END


      SUBROUTINE MATMUL_3X3(A, B, C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
      INTEGER I, J, K
C
      DO I = 1, 3
        DO J = 1, 3
          C(I,J) = 0.0D0
          DO K = 1, 3
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END


C ============================================================
C  KMATINV3 - 3×3 inverse (from mat_utils.f, duplicated for
C  self-contained compilation within umat_user.f)
C ============================================================
      SUBROUTINE KMATINV3(A, AINV)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), AINV(3,3)
      DOUBLE PRECISION DET, EM20
      PARAMETER (EM20=1.0D-20)
C
      DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1    - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2    + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      IF(ABS(DET) .LT. EM20) DET = EM20
C
      AINV(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
      AINV(2,1) = (A(3,2)*A(1,3)-A(1,2)*A(3,3))/DET
      AINV(3,1) = (A(1,2)*A(2,3)-A(2,2)*A(1,3))/DET
      AINV(1,2) = (A(2,3)*A(3,1)-A(2,1)*A(3,3))/DET
      AINV(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
      AINV(3,2) = (A(2,1)*A(1,3)-A(1,1)*A(2,3))/DET
      AINV(1,3) = (A(2,1)*A(3,2)-A(3,1)*A(2,2))/DET
      AINV(2,3) = (A(3,1)*A(1,2)-A(1,1)*A(3,2))/DET
      AINV(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
C
      RETURN
      END


C ============================================================
C  PRODMAT - 3×3 matrix multiply C = A * B
C ============================================================
      SUBROUTINE PRODMAT(A, B, C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
      INTEGER I, J, K
C
      DO I = 1, 3
        DO J = 1, 3
          C(I,J) = 0.0D0
          DO K = 1, 3
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END


C ============================================================
C  PRODAAT - symmetric product C = A * A^T
C ============================================================
      SUBROUTINE PRODAAT(A, C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), C(3,3)
C
      C(1,1) = A(1,1)*A(1,1)+A(1,2)*A(1,2)+A(1,3)*A(1,3)
      C(1,2) = A(1,1)*A(2,1)+A(1,2)*A(2,2)+A(1,3)*A(2,3)
      C(1,3) = A(1,1)*A(3,1)+A(1,2)*A(3,2)+A(1,3)*A(3,3)
      C(2,1) = C(1,2)
      C(2,2) = A(2,1)*A(2,1)+A(2,2)*A(2,2)+A(2,3)*A(2,3)
      C(2,3) = A(2,1)*A(3,1)+A(2,2)*A(3,2)+A(2,3)*A(3,3)
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)
      C(3,3) = A(3,1)*A(3,1)+A(3,2)*A(3,2)+A(3,3)*A(3,3)
C
      RETURN
      END
