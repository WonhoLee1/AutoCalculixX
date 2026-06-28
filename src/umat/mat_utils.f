C ============================================================
C  WHTOOLs MATCALIB 2026 — mat_utils.f
C  3×3 matrix operations for PRF UMAT
C  Ported from OpenRadioss:
C    kmatinv3.F  - 3×3 matrix inverse
C    prodmat.F   - C = A * B (3×3 multiply)
C    prodaat.F   - C = A * A^T (symmetric)
C    calcmatb    - combined: fe = f * inv(fp), matb = fe * fe^T
C ============================================================

C ============================================================
C  KMATINV3 - Inverse of a 3×3 matrix
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
      AINV(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / DET
      AINV(2,1) = (A(3,2)*A(1,3) - A(1,2)*A(3,3)) / DET
      AINV(3,1) = (A(1,2)*A(2,3) - A(2,2)*A(1,3)) / DET
      AINV(1,2) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / DET
      AINV(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / DET
      AINV(3,2) = (A(2,1)*A(1,3) - A(1,1)*A(2,3)) / DET
      AINV(1,3) = (A(2,1)*A(3,2) - A(3,1)*A(2,2)) / DET
      AINV(2,3) = (A(3,1)*A(1,2) - A(1,1)*A(3,2)) / DET
      AINV(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / DET
C
      RETURN
      END


C ============================================================
C  PRODMAT - 3×3 matrix multiply: C = A * B
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
C  PRODAAT - Symmetric product: C = A * A^T
C  (left Cauchy-Green tensor from deformation gradient)
C ============================================================
      SUBROUTINE PRODAAT(A, C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), C(3,3)
C
      C(1,1) = A(1,1)*A(1,1) + A(1,2)*A(1,2) + A(1,3)*A(1,3)
      C(1,2) = A(1,1)*A(2,1) + A(1,2)*A(2,2) + A(1,3)*A(2,3)
      C(1,3) = A(1,1)*A(3,1) + A(1,2)*A(3,2) + A(1,3)*A(3,3)
      C(2,1) = C(1,2)
      C(2,2) = A(2,1)*A(2,1) + A(2,2)*A(2,2) + A(2,3)*A(2,3)
      C(2,3) = A(2,1)*A(3,1) + A(2,2)*A(3,2) + A(2,3)*A(3,3)
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)
      C(3,3) = A(3,1)*A(3,1) + A(3,2)*A(3,2) + A(3,3)*A(3,3)
C
      RETURN
      END


C ============================================================
C  CALCMATB - Compute elastic left Cauchy-Green:
C    fe = f * inv(fp)
C    matb = fe * fe^T
C
C  This is the trial elastic "b" for a viscous network.
C ============================================================
      SUBROUTINE CALCMATB(F, FP, MATB)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), FP(3,3), MATB(3,3)
      DOUBLE PRECISION INVFP(3,3), FE(3,3)
C
      CALL KMATINV3(FP, INVFP)
      CALL PRODMAT(F, INVFP, FE)
      CALL PRODAAT(FE, MATB)
C
      RETURN
      END
