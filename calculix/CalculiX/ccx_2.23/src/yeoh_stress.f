C ============================================================
C  WHTOOLs MATCALIB 2026 — yeoh_stress.f
C  Yeoh (Reduced Polynomial 3rd order) hyperelastic stress
C
C  Strain energy (deviatoric + volumetric split):
C    W_dev  = C10*(I1b-3) + C20*(I1b-3)^2 + C30*(I1b-3)^3
C    U(J)   = D1*(J-1)^2 + D2*(J-1)^4 + D3*(J-1)^6
C
C  where I1b = J^(-2/3) * I1 (deviatoric first invariant)
C
C  Input:  MATB(3,3) = left Cauchy-Green b = F*F^T
C  Output: SIG(3,3)  = Cauchy stress
C ============================================================
      SUBROUTINE YEOH_STRESS(
     1     MATB, C10, C20, C30, D1, D2, D3,
     2     SIG, BI1, JDET)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION MATB(3,3), C10, C20, C30, D1, D2, D3
      DOUBLE PRECISION SIG(3,3), BI1, JDET
C
      DOUBLE PRECISION I1, I1B, JTHIRD, J2THIRD, I1B_M3
      DOUBLE PRECISION DPHIDI1, DPHIDJ, DJ, DJ2, DJ3
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, EM20, THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, SIX=6.0D0, EM20=1.0D-20)
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
      I1B_M3 = I1B - THREE
      DPHIDI1 = C10
     1        + TWO  * C20 * I1B_M3
     2        + THREE * C30 * I1B_M3 * I1B_M3
      DPHIDI1 = TWO * DPHIDI1 / MAX(EM20, JDET)
C     (includes 2/J factor like Arruda-Boyce)
C
C --- dU/dJ (volumetric) ---
C     U = D1*(J-1)^2 + D2*(J-1)^4 + D3*(J-1)^6
C     dU/dJ = 2*D1*(J-1) + 4*D2*(J-1)^3 + 6*D3*(J-1)^5
      DJ  = JDET - ONE
      DJ2 = DJ * DJ
      DJ3 = DJ2 * DJ
      DPHIDJ = TWO*D1*DJ + FOUR*D2*DJ3 + SIX*D3*DJ3*DJ2
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
