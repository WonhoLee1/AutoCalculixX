C ============================================================
C  WHTOOLs MATCALIB 2026 — polynomial_stress.f
C  Polynomial hyperelastic stress from left Cauchy-Green tensor
C  Ported from OpenRadioss sigpoly.F (POLYSTRESS2)
C
C  Strain energy:
C    W = sum_{i+j=1..3} Cij * (I1_bar - 3)^i * (I2_bar - 3)^j
C      + sum_{k=1..3} Dk * (J - 1)^{2k}
C
C  Special cases:
C    C10  only        => Neo-Hookean
C    C10 + C01        => Mooney-Rivlin (2-parameter)
C    All Cij nonzero  => Full Polynomial
C ============================================================
      SUBROUTINE POLY_STRESS(
     1     MATB, C10, C01, C20, C11, C02, C30, C21, C12, C03,
     2     D1, D2, D3, SIG, BI1, BI2, AJ, IHYPER, RBULK)
C
C  INPUT:
C    MATB(3,3)  : left Cauchy-Green tensor b = F * F^T
C    C10..C03   : deviatoric polynomial coefficients
C    D1..D3     : volumetric coefficients
C    IHYPER     : 1 = D1/D2/D3 form, 2 = bulk modulus form
C    RBULK      : bulk modulus (used when IHYPER=2)
C
C  OUTPUT:
C    SIG(3,3)   : Cauchy stress
C    BI1        : deviatoric first invariant I1_bar
C    BI2        : deviatoric second invariant I2_bar
C    AJ         : volume ratio J = det(F)
C ============================================================
      IMPLICIT NONE

      DOUBLE PRECISION MATB(3,3), C10, C01, C20, C11, C02
      DOUBLE PRECISION C30, C21, C12, C03, D1, D2, D3
      DOUBLE PRECISION SIG(3,3), BI1, BI2, AJ, RBULK
      INTEGER IHYPER

      DOUBLE PRECISION MATB2(3,3)
      DOUBLE PRECISION I1, I2, TRB2, TRB22
      DOUBLE PRECISION JTHIRD, J2THIRD, J4THIRD
      DOUBLE PRECISION DPHIDI1, DPHIDI2, DPHIDJ
      DOUBLE PRECISION AA, BB, CC, INV2J
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, EM20
      DOUBLE PRECISION THIRD, TWO_THIRD, FOUR_THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, SIX=6.0D0, EM20=1.0D-20)
      PARAMETER (THIRD=1.0D0/3.0D0, TWO_THIRD=2.0D0/3.0D0,
     1           FOUR_THIRD=4.0D0/3.0D0)

C --- compute MATB2 = MATB * MATB ---
      MATB2(1,1) = MATB(1,1)*MATB(1,1)+MATB(1,2)*MATB(2,1)
     1           + MATB(1,3)*MATB(3,1)
      MATB2(1,2) = MATB(1,1)*MATB(1,2)+MATB(1,2)*MATB(2,2)
     1           + MATB(1,3)*MATB(3,2)
      MATB2(1,3) = MATB(1,1)*MATB(1,3)+MATB(1,2)*MATB(2,3)
     1           + MATB(1,3)*MATB(3,3)
      MATB2(2,1) = MATB2(1,2)
      MATB2(2,2) = MATB(2,1)*MATB(1,2)+MATB(2,2)*MATB(2,2)
     1           + MATB(2,3)*MATB(3,2)
      MATB2(2,3) = MATB(2,1)*MATB(1,3)+MATB(2,2)*MATB(2,3)
     1           + MATB(2,3)*MATB(3,3)
      MATB2(3,1) = MATB2(1,3)
      MATB2(3,2) = MATB2(2,3)
      MATB2(3,3) = MATB(3,1)*MATB(1,3)+MATB(3,2)*MATB(2,3)
     1           + MATB(3,3)*MATB(3,3)

C --- J = sqrt(det(b)) ---
      AJ = MATB(1,1)*(MATB(2,2)*MATB(3,3)-MATB(2,3)*MATB(3,2))
     1   - MATB(1,2)*(MATB(2,1)*MATB(3,3)-MATB(2,3)*MATB(3,1))
     2   + MATB(1,3)*(MATB(2,1)*MATB(3,2)-MATB(2,2)*MATB(3,1))
      AJ = SQRT(MAX(EM20, AJ))

C --- invariants ---
      I1 = MATB(1,1) + MATB(2,2) + MATB(3,3)
      TRB2  = I1 * I1
      TRB22 = MATB2(1,1) + MATB2(2,2) + MATB2(3,3)
      I2 = (TRB2 - TRB22) / TWO

C --- deviatoric invariants ---
      IF (AJ .GT. ZERO) THEN
        JTHIRD  = EXP(-THIRD * LOG(AJ))
        J2THIRD = JTHIRD * JTHIRD
        J4THIRD = J2THIRD * J2THIRD
      ELSE
        JTHIRD  = ZERO
        J2THIRD = ZERO
        J4THIRD = ZERO
      END IF
      BI1 = I1 * J2THIRD
      BI2 = I2 * J4THIRD

C --- derivatives of strain energy ---
      DPHIDI1 = C10
     1        + TWO*C20*(BI1-THREE)
     2        + THREE*C30*(BI1-THREE)*(BI1-THREE)
     3        + C11*(BI2-THREE)
     4        + C12*(BI2-THREE)*(BI2-THREE)
     5        + TWO*C21*(BI1-THREE)*(BI2-THREE)

      DPHIDI2 = C01
     1        + TWO*C02*(BI2-THREE)
     2        + THREE*C03*(BI2-THREE)*(BI2-THREE)
     3        + C11*(BI1-THREE)
     4        + C21*(BI1-THREE)*(BI1-THREE)
     5        + TWO*C12*(BI1-THREE)*(BI2-THREE)

C --- volumetric pressure ---
      IF (IHYPER .EQ. 1) THEN
        DPHIDJ = TWO*D1*(AJ-ONE)
     1         + FOUR*D2*(AJ-ONE)*(AJ-ONE)*(AJ-ONE)
     2         + SIX*D3*(AJ-ONE)*(AJ-ONE)*(AJ-ONE)
     3              *(AJ-ONE)*(AJ-ONE)
      ELSE
        DPHIDJ = RBULK * (ONE - ONE/MAX(EM20,AJ))
      END IF

      INV2J = TWO / MAX(EM20, AJ)

C --- Cauchy stress: sigma = (2/J)[(dW/dI1 + dW/dI2*I1)*b - dW/dI2*b^2]
C     - (2/3J)(I1*dW/dI1 + 2*I2*dW/dI2)*I + dW/dJ*I ---
      AA = (DPHIDI1 + DPHIDI2*BI1) * INV2J * J2THIRD
      BB = DPHIDI2 * INV2J * J4THIRD
      CC = THIRD * INV2J * (BI1*DPHIDI1 + TWO*BI2*DPHIDI2)

      SIG(1,1) = AA*MATB(1,1) - BB*MATB2(1,1) - CC + DPHIDJ
      SIG(2,2) = AA*MATB(2,2) - BB*MATB2(2,2) - CC + DPHIDJ
      SIG(3,3) = AA*MATB(3,3) - BB*MATB2(3,3) - CC + DPHIDJ
      SIG(1,2) = AA*MATB(1,2) - BB*MATB2(1,2)
      SIG(1,3) = AA*MATB(1,3) - BB*MATB2(1,3)
      SIG(2,3) = AA*MATB(2,3) - BB*MATB2(2,3)
      SIG(2,1) = SIG(1,2)
      SIG(3,1) = SIG(1,3)
      SIG(3,2) = SIG(2,3)

      RETURN
      END

C ============================================================
C  Compute left Cauchy-Green tensor b from deformation gradient F
C ============================================================
      SUBROUTINE CALC_MATB(F, B)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), B(3,3)

      B(1,1) = F(1,1)*F(1,1) + F(1,2)*F(1,2) + F(1,3)*F(1,3)
      B(1,2) = F(1,1)*F(2,1) + F(1,2)*F(2,2) + F(1,3)*F(2,3)
      B(1,3) = F(1,1)*F(3,1) + F(1,2)*F(3,2) + F(1,3)*F(3,3)
      B(2,1) = B(1,2)
      B(2,2) = F(2,1)*F(2,1) + F(2,2)*F(2,2) + F(2,3)*F(2,3)
      B(2,3) = F(2,1)*F(3,1) + F(2,2)*F(3,2) + F(2,3)*F(3,3)
      B(3,1) = B(1,3)
      B(3,2) = B(2,3)
      B(3,3) = F(3,1)*F(3,1) + F(3,2)*F(3,2) + F(3,3)*F(3,3)

      RETURN
      END
