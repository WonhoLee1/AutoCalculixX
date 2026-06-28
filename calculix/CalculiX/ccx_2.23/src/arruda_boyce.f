C ============================================================
C  arruda_boyce.f
C  Arruda-Boyce 8-chain hyperelastic Cauchy stress
C  Ported from OpenRadioss sigaboyce.F
C
C  Strain energy:
C    W = MU * sum_{i=1..5} Ci * beta^(i-1) * (I1_bar^i - 3^i)
C      + (1/2)*D * (J - 1)^2
C
C  where:
C    C1 = 1/2, C2 = 1/20, C3 = 11/1050, C4 = 19/7000, C5 = 519/673750
C    (Arruda-Boyce 8-chain coefficients)
C
C  Input: MATB = left Cauchy-Green b = F * F^T
C  Output: SIG(3,3) = Cauchy stress
C ============================================================
      SUBROUTINE ARRUDA_BOYCE(
     1     MATB, C1, C2, C3, C4, C5, MU, D, BETA,
     2     SIG, BI1, JDET)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION MATB(3,3), C1, C2, C3, C4, C5, MU, D, BETA
      DOUBLE PRECISION SIG(3,3), BI1, JDET
C
      DOUBLE PRECISION I1
      DOUBLE PRECISION JTHIRD, J2THIRD
      DOUBLE PRECISION DPHIDI1, DPHIDJ
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, FIVE, EM20, THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, FIVE=5.0D0, EM20=1.0D-20)
      PARAMETER (THIRD=1.0D0/3.0D0)
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
C --- deviatoric first invariant I1_bar ---
      IF (JDET .GT. ZERO) THEN
        JTHIRD  = EXP((-THIRD) * LOG(JDET))
        J2THIRD = JTHIRD * JTHIRD
      ELSE
        JTHIRD  = ZERO
        J2THIRD = ZERO
      END IF
      BI1 = I1 * J2THIRD
C
C --- dW/dI1_bar (includes 2/J factor like OpenRadioss) ---
C     Note: DPHIDI1 = (2/J) * dW/dI1_bar
      DPHIDI1 = TWO * MU * ( C1
     1     + TWO   * C2 * BETA          * BI1
     2     + THREE * C3 * (BETA*BI1)**2
     3     + FOUR  * C4 * (BETA*BI1)**3
     4     + FIVE  * C5 * (BETA*BI1)**4 ) / MAX(EM20, JDET)
C
C --- dW/dJ (volumetric) ---
C     D in the caller is 1/D (i.e. D = 1/D1 in polynomial convention)
      DPHIDJ = D * (JDET - ONE / MAX(EM20, JDET))
C
C --- Cauchy stress (same formulation as OpenRadioss sigaboyce) ---
C     DIAGONAL: sigma_ii = dW/dI1*(b_ii - 1/3*I1_bar) + dW/dJ
      SIG(1,1) = DPHIDI1 * (MATB(1,1) - THIRD*BI1) + DPHIDJ
      SIG(2,2) = DPHIDI1 * (MATB(2,2) - THIRD*BI1) + DPHIDJ
      SIG(3,3) = DPHIDI1 * (MATB(3,3) - THIRD*BI1) + DPHIDJ
C     OFF-DIAGONAL: sigma_ij = dW/dI1 * b_ij
      SIG(1,2) = DPHIDI1 * MATB(1,2)
      SIG(2,3) = DPHIDI1 * MATB(2,3)
      SIG(3,1) = DPHIDI1 * MATB(3,1)
      SIG(2,1) = SIG(1,2)
      SIG(3,2) = SIG(2,3)
      SIG(1,3) = SIG(3,1)
C
      RETURN
      END
