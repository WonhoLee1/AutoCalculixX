C ============================================================
C  WHTOOLs MATCALIB 2026 — ogden_stress.f
C  Ogden hyperelastic stress (N=3 terms)
C
C  Strain energy:
C    W_dev  = sum_{k=1..3} (2*mu_k/alpha_k^2) *
C             (lam1b^alpha_k + lam2b^alpha_k + lam3b^alpha_k - 3)
C    U(J)   = D1*(J-1)^2 + D2*(J-1)^4 + D3*(J-1)^6
C
C  where lamb_i = deviatoric principal stretches = J^(-1/3) * sqrt(eig_i(b))
C
C  Input:  MATB(3,3) = left Cauchy-Green b = F*F^T
C  Output: SIG(3,3)  = Cauchy stress
C ============================================================
      SUBROUTINE OGDEN_STRESS(
     1     MATB, MU1, ALP1, MU2, ALP2, MU3, ALP3,
     2     D1, D2, D3, SIG, BI1, JDET)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION MATB(3,3)
      DOUBLE PRECISION MU1, ALP1, MU2, ALP2, MU3, ALP3
      DOUBLE PRECISION D1, D2, D3
      DOUBLE PRECISION SIG(3,3), BI1, JDET
C
      DOUBLE PRECISION B(3,3), EVEC(3,3), EVAL(3)
      DOUBLE PRECISION I1, JTHIRD, LB(3), L2(3)
      DOUBLE PRECISION DJ, DJ2, DJ3, DPHIDJ
      DOUBLE PRECISION SDEV(3)
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, SIX, EM20
      DOUBLE PRECISION THIRD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, SIX=6.0D0, EM20=1.0D-20)
      PARAMETER (THIRD=1.0D0/3.0D0)
C
      INTEGER I, J, K
C
C --- J = sqrt(det(b)) ---
      JDET = MATB(1,1)*(MATB(2,2)*MATB(3,3)-MATB(2,3)*MATB(3,2))
     1     - MATB(1,2)*(MATB(2,1)*MATB(3,3)-MATB(2,3)*MATB(3,1))
     2     + MATB(1,3)*(MATB(2,1)*MATB(3,2)-MATB(2,2)*MATB(3,1))
      JDET = SQRT(MAX(EM20, JDET))
C
C --- I1 = trace(b) ---
      I1 = MATB(1,1) + MATB(2,2) + MATB(3,3)
      BI1 = I1 * EXP((-2.0D0/3.0D0)*LOG(MAX(EM20,JDET)))
C
C --- copy b to local symmetric array ---
      DO I = 1, 3
        DO J = 1, 3
          B(I,J) = MATB(I,J)
        END DO
      END DO
C
C --- eigenvalue decomposition of b (symmetric) via Jacobi ---
      CALL JACOBI3(B, EVEC, EVAL, 3, 50)
C
C --- eigenvalues (principal stretches squared) sorted EVAL(1) >= EVAL(2) >= EVAL(3) ---
C     lambda_i = sqrt(eval_i), then deviatoric: lb_i = J^(-1/3) * lambda_i
      IF (JDET .GT. ZERO) THEN
        JTHIRD = EXP((-THIRD) * LOG(JDET))
      ELSE
        JTHIRD = ZERO
      END IF
      DO K = 1, 3
        L2(K) = MAX(EM20, EVAL(K))
        LB(K) = SQRT(L2(K)) * JTHIRD
      END DO
C
C --- principal deviatoric Cauchy stress (without pressure) ---
C     S_i = (1/J) * sum_k (2*mu_k/alpha_k) * (lb_i^alpha_k)
C     minus volumetric pressure dU/dJ
C
      SDEV(1) = ZERO
      SDEV(2) = ZERO
      SDEV(3) = ZERO
C
      IF (ABS(ALP1) .GT. EM20) THEN
        SDEV(1) = SDEV(1) + (TWO*MU1/ALP1) * LB(1)**ALP1
        SDEV(2) = SDEV(2) + (TWO*MU1/ALP1) * LB(2)**ALP1
        SDEV(3) = SDEV(3) + (TWO*MU1/ALP1) * LB(3)**ALP1
      ELSE
C       alpha -> 0: limit is MU*ln(lambda)
        SDEV(1) = SDEV(1) + TWO*MU1 * LOG(MAX(EM20, LB(1)))
        SDEV(2) = SDEV(2) + TWO*MU1 * LOG(MAX(EM20, LB(2)))
        SDEV(3) = SDEV(3) + TWO*MU1 * LOG(MAX(EM20, LB(3)))
      END IF
C
      IF (ABS(ALP2) .GT. EM20) THEN
        SDEV(1) = SDEV(1) + (TWO*MU2/ALP2) * LB(1)**ALP2
        SDEV(2) = SDEV(2) + (TWO*MU2/ALP2) * LB(2)**ALP2
        SDEV(3) = SDEV(3) + (TWO*MU2/ALP2) * LB(3)**ALP2
      ELSE
        SDEV(1) = SDEV(1) + TWO*MU2 * LOG(MAX(EM20, LB(1)))
        SDEV(2) = SDEV(2) + TWO*MU2 * LOG(MAX(EM20, LB(2)))
        SDEV(3) = SDEV(3) + TWO*MU2 * LOG(MAX(EM20, LB(3)))
      END IF
C
      IF (ABS(ALP3) .GT. EM20) THEN
        SDEV(1) = SDEV(1) + (TWO*MU3/ALP3) * LB(1)**ALP3
        SDEV(2) = SDEV(2) + (TWO*MU3/ALP3) * LB(2)**ALP3
        SDEV(3) = SDEV(3) + (TWO*MU3/ALP3) * LB(3)**ALP3
      ELSE
        SDEV(1) = SDEV(1) + TWO*MU3 * LOG(MAX(EM20, LB(1)))
        SDEV(2) = SDEV(2) + TWO*MU3 * LOG(MAX(EM20, LB(2)))
        SDEV(3) = SDEV(3) + TWO*MU3 * LOG(MAX(EM20, LB(3)))
      END IF
C
C --- principal deviatoric Cauchy stress (unscaled, before hydrostatic subtraction) ---
C     S_full_i = (1/J) * sum_k (2*mu_k/alpha_k) * lb_i^alpha_k
      DO K = 1, 3
        SDEV(K) = SDEV(K) / MAX(EM20, JDET)
      END DO
C
C --- subtract hydrostatic part to make deviatoric ---
C     S_dev_i = S_full_i - (1/3)*(S_full_1 + S_full_2 + S_full_3)
      DPHIDJ = THIRD * (SDEV(1) + SDEV(2) + SDEV(3))
      SDEV(1) = SDEV(1) - DPHIDJ
      SDEV(2) = SDEV(2) - DPHIDJ
      SDEV(3) = SDEV(3) - DPHIDJ
C
C --- volumetric pressure: dU/dJ (added back to deviatoric) ---
C     U = D1*(J-1)^2 + D2*(J-1)^4 + D3*(J-1)^6
C     dU/dJ = 2*D1*(J-1) + 4*D2*(J-1)^3 + 6*D3*(J-1)^5
      DJ  = JDET - ONE
      DJ2 = DJ * DJ
      DJ3 = DJ2 * DJ
      DPHIDJ = TWO*D1*DJ + FOUR*D2*DJ3 + SIX*D3*DJ3*DJ2
C
C --- total principal stress = deviatoric + volumetric ---
      SDEV(1) = SDEV(1) + DPHIDJ
      SDEV(2) = SDEV(2) + DPHIDJ
      SDEV(3) = SDEV(3) + DPHIDJ
C
C --- rotate back to Cartesian: SIG = V * diag(SDEV) * V^T ---
      DO I = 1, 3
        DO J = 1, 3
          SIG(I,J) = ZERO
        END DO
      END DO
      DO K = 1, 3
        DO I = 1, 3
          DO J = 1, 3
            SIG(I,J) = SIG(I,J) + EVEC(I,K) * SDEV(K) * EVEC(J,K)
          END DO
        END DO
      END DO
C
C --- enforce symmetry ---
      SIG(2,1) = SIG(1,2)
      SIG(3,2) = SIG(2,3)
      SIG(1,3) = SIG(3,1)
C
      RETURN
      END
C
C
C ============================================================
C  JACOBI3 - Symmetric Jacobi eigenvalue solver for 3x3 matrix
C  Self-contained, fully updates both triangles for symmetry.
C
C  Input:  A(3,3)    = symmetric input matrix
C          N         = matrix dimension (must be 3)
C          MAX_SWEEP = maximum Jacobi sweeps
C  Output: VEC(3,3)  = eigenvectors (columns, VEC(:,i) is eigenvector i)
C          EVAL(3)   = eigenvalues (sorted descending), EVAL(1) >= EVAL(2) >= EVAL(3)
C ============================================================
      SUBROUTINE JACOBI3(A, VEC, EVAL, N, MAX_SWEEP)
      IMPLICIT NONE
      INTEGER N, MAX_SWEEP
      DOUBLE PRECISION A(3,3), VEC(3,3), EVAL(3)
C
      DOUBLE PRECISION B(3,3), V(3,3)
      DOUBLE PRECISION THETA, T, C, S, TAU
      DOUBLE PRECISION TMP, AKK, AII, AIKAK
      DOUBLE PRECISION SM, SMIN, DTOL
      INTEGER I, J, K, IP, IQ, SWEEP
      DOUBLE PRECISION ZERO, ONE, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, HALF=0.5D0)
C
C --- initialise B = A (symmetric), V = identity ---
      DO I = 1, 3
        DO J = 1, 3
          B(I,J) = A(I,J)
          V(I,J) = ZERO
        END DO
        V(I,I) = ONE
      END DO
C
      SMIN = 1.0D-15
C
      DO SWEEP = 1, MAX_SWEEP
C       --- compute sum of off-diagonals (half the matrix, since symmetric) ---
        SM = ZERO
        DO IP = 1, N - 1
          DO IQ = IP + 1, N
            SM = SM + ABS(B(IP,IQ))
          END DO
        END DO
        IF (SM .LT. SMIN) EXIT
C
        DO IP = 1, N - 1
          DO IQ = IP + 1, N
            IF (ABS(B(IP,IQ)) .LT. SMIN) CYCLE
C
C           --- compute Jacobi rotation angle ---
            IF (ABS(B(IQ,IQ) - B(IP,IP)) .LT. SMIN) THEN
              THETA = HALF * SIGN(ONE, B(IP,IQ))
            ELSE
              THETA = HALF * (B(IQ,IQ) - B(IP,IP)) / B(IP,IQ)
            END IF
            T = SIGN(ONE /(ABS(THETA) + SQRT(ONE + THETA*THETA)), THETA)
            C = ONE / SQRT(ONE + T*T)
            S = T * C
            TAU = S / (ONE + C)
C
C           --- update diagonal of B ---
            B(IP,IP) = B(IP,IP) - T * B(IP,IQ)
            B(IQ,IQ) = B(IQ,IQ) + T * B(IP,IQ)
            B(IP,IQ) = ZERO
            B(IQ,IP) = ZERO
C
C           --- update remaining rows/columns (preserve symmetry) ---
            DO J = 1, N
              IF (J .NE. IP .AND. J .NE. IQ) THEN
                AKK = B(J,IP)
                AII = B(J,IQ)
                B(J,IP) = C*AKK - S*AII
                B(IP,J) = B(J,IP)
                B(J,IQ) = S*AKK + C*AII
                B(IQ,J) = B(J,IQ)
              END IF
            END DO
C
C           --- update eigenvector matrix ---
            DO J = 1, N
              AKK = V(J,IP)
              AII = V(J,IQ)
              V(J,IP) = C*AKK - S*AII
              V(J,IQ) = S*AKK + C*AII
            END DO
          END DO
        END DO
      END DO
C
C --- copy diagonal to eigenvalues ---
      DO I = 1, N
        EVAL(I) = B(I,I)
      END DO
C
C --- sort eigenvalues descending with corresponding eigenvector columns ---
      DO I = 1, N - 1
        DO J = I + 1, N
          IF (EVAL(J) .GT. EVAL(I)) THEN
            TMP = EVAL(I)
            EVAL(I) = EVAL(J)
            EVAL(J) = TMP
            DO K = 1, N
              TMP = V(K,I)
              V(K,I) = V(K,J)
              V(K,J) = TMP
            END DO
          END IF
        END DO
      END DO
C
C --- copy V to VEC output ---
      DO I = 1, 3
        DO J = 1, 3
          VEC(I,J) = V(I,J)
        END DO
      END DO
C
      RETURN
      END