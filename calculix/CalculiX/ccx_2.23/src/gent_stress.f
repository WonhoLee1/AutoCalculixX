C ============================================================
C  WHTOOLs MATCALIB 2026 — gent_stress.f
C  Gent hyperelastic stress (locking model)
C
C  Strain energy:
C    W_dev  = -(MU/2) * Jm * ln(1 - (I1b-3)/Jm)
C    U(J)   = D * (J-1)^2   (via dU/dJ = D*(J - 1/J))
C
C  where I1b = J^(-2/3) * I1 (deviatoric first invariant)
C        Jm  = locking parameter (max I1b-3)
C
C  Input:  MATB(3,3) = left Cauchy-Green b = F*F^T
C  Output: SIG(3,3)  = Cauchy stress
C ============================================================
      SUBROUTINE GENT_STRESS(
     1     MATB, MU, JM, D,
     2     SIG, BI1, JDET)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION MATB(3,3), MU, JM, D
      DOUBLE PRECISION SIG(3,3), BI1, JDET
C
      DOUBLE PRECISION I1, I1B, JTHIRD, J2THIRD, ARG
      DOUBLE PRECISION DPHIDI1, DPHIDJ
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, EM20, THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (EM20=1.0D-20, THIRD=1.0D0/3.0D0)
C
C --- J = sqrt(det(b)) ---
      JDET = MATB(1,1)*(MATB(2,2)*MATB(3,3)-MATB(2,3)*MATB(3,2))
     1     - MATB(1,2)*(MATB(2,1)*MATB(3,3)-MATB(2,3)*MATB(3,1))
     2     + MATB(1,3)*(MATB(2,1)*MATB(3,2)-MATB(2,2)*MATB(3,1))
      JDET = SQRT(MAX(EM20, JDET))
C
C --- I1 = trace(b) ---
      I1 = MATB(1,1) + MATB(2,2) + MATB(3,3)
C
C --- deviatoric first invariant I1b = J^(-2/3) * I1 ---
      IF (JDET .GT. ZERO) THEN
        JTHIRD  = EXP((-THIRD) * LOG(JDET))
        J2THIRD = JTHIRD * JTHIRD
      ELSE
        JTHIRD  = ZERO
        J2THIRD = ZERO
      END IF
      I1B = I1 * J2THIRD
      BI1 = I1B
C
C --- dW_dev/dI1b ---
C     W = -(MU/2)*Jm*ln(1 - (I1b-3)/Jm)
C     dW/dI1b = (MU/2) / (1 - (I1b-3)/Jm) = MU*Jm / (2*(Jm - (I1b-3)))
      ARG = JM - (I1B - THREE)
      IF (ARG .LE. EM20) THEN
C       Locking reached — use large stiffness as penalty
        DPHIDI1 = HUGE(ONE) * 1.0D-10
      ELSE
        DPHIDI1 = MU * JM / (TWO * ARG)
      END IF
      DPHIDI1 = TWO * DPHIDI1 / MAX(EM20, JDET)
C     (includes 2/J factor)
C
C --- dU/dJ (volumetric) ---
C     Same as Arruda-Boyce: DPHIDJ = D*(J - 1/J)
      DPHIDJ = D * (JDET - ONE / MAX(EM20, JDET))
C
C --- Cauchy stress ---
      SIG(1,1) = DPHIDI1 * (MATB(1,1) - THIRD*I1B) + DPHIDJ
      SIG(2,2) = DPHIDI1 * (MATB(2,2) - THIRD*I1B) + DPHIDJ
      SIG(3,3) = DPHIDI1 * (MATB(3,3) - THIRD*I1B) + DPHIDJ
      SIG(1,2) = DPHIDI1 * MATB(1,2)
      SIG(2,3) = DPHIDI1 * MATB(2,3)
      SIG(3,1) = DPHIDI1 * MATB(3,1)
      SIG(2,1) = SIG(1,2)
      SIG(3,2) = SIG(2,3)
      SIG(1,3) = SIG(3,1)
C
      RETURN
      END
