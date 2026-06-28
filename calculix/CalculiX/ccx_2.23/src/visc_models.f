C ============================================================
C  visc_models.f
C  Viscous creep models for PRF (Parallel Rheological Framework)
C  Ported from OpenRadioss:
C    viscbb.F    - Bergstrom-Boyce
C    viscpower.F - Power law
C    viscsinh.F  - Hyperbolic sine
C
C  All subroutines compute the creep strain increment DGAMMA.
C ============================================================

C ============================================================
C  VISC_BB - Bergstrom-Boyce creep model
C
C  dgamma = A1 * (lambda - 1 + KSI)^EXPC * (tau_norm / tauref)^EXPM
C  where  lambda = sqrt(1/3 * (Fp_11^2 + Fp_22^2 + Fp_33^2))
C
C  Input:
C    FP(3,3)   : viscous deformation gradient (trial)
C    TBNORM    : effective stress norm ||dev(sigma_v)||
C    A1        : pre-exponential factor (creep rate coefficient)
C    EXPC      : chain stretch exponent
C    EXPM      : stress exponent
C    KSI       : chain stretch offset
C    TAUREF    : reference stress
C  Output:
C    DGAMMA    : creep strain increment
C ============================================================
      SUBROUTINE VISC_BB(FP, TBNORM, A1, EXPC, EXPM, KSI,
     1                   TAUREF, DGAMMA)
      IMPLICIT NONE
      DOUBLE PRECISION FP(3,3), TBNORM, A1, EXPC, EXPM, KSI
      DOUBLE PRECISION TAUREF, DGAMMA
C
      DOUBLE PRECISION IP1, LPCHAIN, TEMP
      DOUBLE PRECISION ZERO, ONE, THREE, THIRD, EM20
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, THREE=3.0D0)
      PARAMETER (THIRD=ONE/THREE, EM20=1.0D-20)
C
C     chain stretch from diagonal of viscous Fp
      IP1 = FP(1,1)*FP(1,1) + FP(2,2)*FP(2,2) + FP(3,3)*FP(3,3)
      LPCHAIN = SQRT(MAX(ZERO, THIRD*IP1))
C
C     dgamma = A1 * (lambda - 1 + KSI)^EXPC * (tau/tauref)^EXPM
      TEMP = MAX(EM20, LPCHAIN - ONE + KSI)
      DGAMMA = A1 * EXP(EXPC * LOG(TEMP))
     1       * (TBNORM / MAX(EM20, TAUREF))**EXPM
C
      RETURN
      END


C ============================================================
C  VISC_POWER - Power law creep
C
C  Integrated form:
C    dgamma = A1 * [ ((expm+1)*gamma_old)^expm * tau_norm^expn ]^(1/(1+expm))
C
C  Input:
C    TBNORM    : effective stress norm
C    A1        : pre-exponential factor
C    EXPM      : strain hardening exponent
C    EXPN      : stress exponent
C    GAMMAOLD  : cumulative equivalent creep strain (beginning of step)
C  Output:
C    DGAMMA    : creep strain increment
C ============================================================
      SUBROUTINE VISC_POWER(TBNORM, A1, EXPM, EXPN,
     1                      GAMMAOLD, DGAMMA)
      IMPLICIT NONE
      DOUBLE PRECISION TBNORM, A1, EXPM, EXPN
      DOUBLE PRECISION GAMMAOLD, DGAMMA
C
      DOUBLE PRECISION TEMP1, TEMP2, TEMP3
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
      TEMP1 = (EXPM + ONE) * GAMMAOLD
      TEMP2 = EXP(EXPN * LOG(TBNORM))
      TEMP3 = EXP(EXPM * LOG(TEMP1))
      DGAMMA = A1 * EXP((ONE/(ONE+EXPM)) * LOG(TEMP2 * TEMP3))
C
      RETURN
      END


C ============================================================
C  VISC_SINH - Hyperbolic sine creep
C
C  dgamma = A1 * sinh(B0 * tau_norm)^EXPN
C
C  Input:
C    TBNORM    : effective stress norm
C    A1        : pre-exponential factor
C    B0        : sinh argument multiplier
C    EXPN      : exponent
C  Output:
C    DGAMMA    : creep strain increment
C ============================================================
      SUBROUTINE VISC_SINH(TBNORM, A1, B0, EXPN, DGAMMA)
      IMPLICIT NONE
      DOUBLE PRECISION TBNORM, A1, B0, EXPN, DGAMMA
C
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
      TEMP = (SINH(B0 * TBNORM))**EXPN
      DGAMMA = A1 * TEMP
C
      RETURN
      END
