C ============================================================
C  WHTOOLs MATCALIB 2026 — hyperelastic_user.f
C  CalculiX UMAT for polynomial hyperelasticity (Neo-Hookean,
C  Mooney-Rivlin, full Polynomial) with numerical tangent
C
C  Ported from OpenRadioss LAW100
C
C  Material constants (elconloc order):
C    elconloc(1) = C10
C    elconloc(2) = C01
C    elconloc(3) = C20
C    elconloc(4) = C11
C    elconloc(5) = C02
C    elconloc(6) = C30
C    elconloc(7) = C21
C    elconloc(8) = C12
C    elconloc(9) = C03
C    elconloc(10)= D1
C    elconloc(11)= D2
C    elconloc(12)= D3
C    elconloc(13)= IHYPER (1=D1/D2/D3, 2=bulk modulus)
C    elconloc(14)= RBULK (bulk modulus, used only if IHYPER=2)
C
C  Usage in .inp:
C    *USER MATERIAL,CONSTANTS=14
C    C10, C01, C20, C11, C02, C30, C21, C12, C03
C    D1,  D2,  D3, 1.0, 0.0
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
C  local variables
      INTEGER NTERMS,NPOLY,IHYPER,NPT,JJ,I,J,ICOL
      PARAMETER (NPOLY=9)
      DOUBLE PRECISION CC(NPOLY),D1,D2,D3,RBULK
      DOUBLE PRECISION F(3,3),B(3,3),SIG(3,3),PK2(3,3)
      DOUBLE PRECISION BI1,BI2,AJ
      DOUBLE PRECISION VEC(3,3),VEC_T(3,3)
      DOUBLE PRECISION FINV(3,3),FINV_T(3,3)
      DOUBLE PRECISION S_REF(6),S_PERT(6)
      DOUBLE PRECISION FPERT(3,3),DELTAF(3,3)
      DOUBLE PRECISION EPS,DET,IWK
      INTEGER K,P,Q,INFO
      PARAMETER (EPS=1.0D-8)
C
C  Voigt index mapping: 1→11, 2→22, 3→33, 4→12, 5→13, 6→23
      INTEGER VPMAP(2,6)
      DATA VPMAP /1,1, 2,2, 3,3, 1,2, 1,3, 2,3/
C
C  ============================================================
C  extract material constants
C  ============================================================
C  number of constants = -kode - 100
      NTERMS = -kode - 100
C
C  default initialization
      DO I = 1, NPOLY
        CC(I) = 0.0D0
      END DO
      D1     = 0.0D0
      D2     = 0.0D0
      D3     = 0.0D0
      IHYPER = 1
      RBULK  = 0.0D0
C
      IF(NTERMS .GE. 1) CC(1)  = elconloc(1)   ! C10
      IF(NTERMS .GE. 2) CC(2)  = elconloc(2)   ! C01
      IF(NTERMS .GE. 3) CC(3)  = elconloc(3)   ! C20
      IF(NTERMS .GE. 4) CC(4)  = elconloc(4)   ! C11
      IF(NTERMS .GE. 5) CC(5)  = elconloc(5)   ! C02
      IF(NTERMS .GE. 6) CC(6)  = elconloc(6)   ! C30
      IF(NTERMS .GE. 7) CC(7)  = elconloc(7)   ! C21
      IF(NTERMS .GE. 8) CC(8)  = elconloc(8)   ! C12
      IF(NTERMS .GE. 9) CC(9)  = elconloc(9)   ! C03
      IF(NTERMS .GE.10) D1     = elconloc(10)  ! D1
      IF(NTERMS .GE.11) D2     = elconloc(11)  ! D2
      IF(NTERMS .GE.12) D3     = elconloc(12)  ! D3
      IF(NTERMS .GE.13) IHYPER = NINT(elconloc(13))
      IF(NTERMS .GE.14) RBULK  = elconloc(14)
C
C  ============================================================
C  copy deformation gradient
C  ============================================================
      DO I = 1, 3
        DO J = 1, 3
          F(I,J) = xkl(J,I)     ! CalculiX stores F transposed
        END DO
      END DO
C
C  ============================================================
C  compute left Cauchy-Green tensor b = F * F^T
C  ============================================================
      CALL CALC_MATB(F, B)
C
C  ============================================================
C  compute Cauchy stress
C  ============================================================
      CALL POLY_STRESS(
     &     B, CC(1), CC(2), CC(3), CC(4), CC(5),
     &     CC(6), CC(7), CC(8), CC(9),
     &     D1, D2, D3, SIG, BI1, BI2, AJ, IHYPER, RBULK)
C
C  ============================================================
C  convert Cauchy stress to PK2 stress: S = J * F^{-1} * σ * F^{-T}
C  ============================================================
C  F^{-1} and F^{-T}
      CALL INV3X3(F, FINV, DET)
      DO I = 1, 3
        DO J = 1, 3
          FINV_T(I,J) = FINV(J,I)
        END DO
      END DO
C
C  PK2 = AJ * FINV * SIG * FINV_T
      CALL MATMUL_3X3(FINV, SIG, VEC)
      CALL MATMUL_3X3(VEC, FINV_T, PK2)
      DO I = 1, 3
        DO J = 1, 3
          PK2(I,J) = AJ * PK2(I,J)
        END DO
      END DO
C
C  store PK2 stress in Voigt form: [11,22,33,12,13,23]
      stre(1) = PK2(1,1)
      stre(2) = PK2(2,2)
      stre(3) = PK2(3,3)
      stre(4) = PK2(1,2)
      stre(5) = PK2(1,3)
      stre(6) = PK2(2,3)
C
C  ============================================================
C  numerical tangent stiffness
C  ============================================================
      IF(icmd .NE. 3) THEN
C
C  save reference PK2 stress (Voigt)
        S_REF(1) = PK2(1,1)
        S_REF(2) = PK2(2,2)
        S_REF(3) = PK2(3,3)
        S_REF(4) = PK2(1,2)
        S_REF(5) = PK2(1,3)
        S_REF(6) = PK2(2,3)
C
C  zero the stiffness array
        DO I = 1, 21
          stiff(I) = 0.0D0
        END DO
C
C  compute FINV_T (needed for perturbation direction)
        DO I = 1, 3
          DO J = 1, 3
            FINV_T(I,J) = FINV(J,I)
          END DO
        END DO
C
C  loop over 6 columns of the tangent
        DO ICOL = 1, 6
          P = VPMAP(1, ICOL)
          Q = VPMAP(2, ICOL)
C
C  construct perturbation direction: H = δ * F^{-T} * (e_p ⊗ e_q)
C  For off-diagonal (p≠q), multiply by 2 to get δE_{pq} = EPS
          DO I = 1, 3
            DO J = 1, 3
              DELTAF(I,J) = 0.0D0
            END DO
          END DO
          IF(P .EQ. Q) THEN
C           normal strain: ΔE_{pp} = EPS
            DO I = 1, 3
              DELTAF(I,P) = EPS * FINV_T(I,P)
            END DO
          ELSE
C           shear strain: need ΔE_{pq} = EPS, but F^{-T}*(e_p⊗e_q) gives EPS/2
C           so multiply by 2
            DO I = 1, 3
              DELTAF(I,P) = 2.0D0 * EPS * FINV_T(I,Q)
            END DO
          END IF
C
C  perturbed deformation gradient: F' = F + δF
          DO I = 1, 3
            DO J = 1, 3
              FPERT(I,J) = F(I,J) + DELTAF(I,J)
            END DO
          END DO
C
C  compute PK2 for perturbed state
          CALL PK2_FROM_F(FPERT, PK2,
     &                    CC(1), CC(2), CC(3), CC(4), CC(5),
     &                    CC(6), CC(7), CC(8), CC(9),
     &                    D1, D2, D3, IHYPER, RBULK)
C
C  finite difference column
          S_PERT(1) = PK2(1,1)
          S_PERT(2) = PK2(2,2)
          S_PERT(3) = PK2(3,3)
          S_PERT(4) = PK2(1,2)
          S_PERT(5) = PK2(1,3)
          S_PERT(6) = PK2(2,3)
C
C  fill Voigt matrix ddsdde(6,6) then convert to stiff(21)
          DO I = 1, 6
            IWK = (S_PERT(I) - S_REF(I)) / EPS
C  store in appropriate stiff position
C  stiff(21) is upper triangle: (1,1),(1,2),(2,2),(1,3),(2,3),(3,3),
C                                 (1,4),(2,4),(3,4),(4,4),
C                                 (1,5),(2,5),(3,5),(4,5),(5,5),
C                                 (1,6),(2,6),(3,6),(4,6),(5,6),(6,6)
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
C  symmetrize off-diagonal terms (average of ∂Si/∂Ej and ∂Sj/∂Ei)
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
C  convergence flag
      pnewdt = -1.0D0
C
      RETURN
      END


C ============================================================
C  PK2_FROM_F: compute PK2 stress from deformation gradient
C  (helper for numerical tangent)
C ============================================================
      SUBROUTINE PK2_FROM_F(F, PK2,
     &     C10, C01, C20, C11, C02, C30, C21, C12, C03,
     &     D1, D2, D3, IHYPER, RBULK)
      IMPLICIT NONE
      DOUBLE PRECISION F(3,3), PK2(3,3)
      DOUBLE PRECISION C10,C01,C20,C11,C02,C30,C21,C12,C03
      DOUBLE PRECISION D1,D2,D3,RBULK
      INTEGER IHYPER
C
      DOUBLE PRECISION B(3,3), SIG(3,3)
      DOUBLE PRECISION FINV(3,3), FINV_T(3,3), VEC(3,3)
      DOUBLE PRECISION BI1, BI2, AJ, DET
      INTEGER I, J
C
      CALL CALC_MATB(F, B)
      CALL POLY_STRESS(B, C10, C01, C20, C11, C02,
     &     C30, C21, C12, C03,
     &     D1, D2, D3, SIG, BI1, BI2, AJ, IHYPER, RBULK)
C
      CALL INV3X3(F, FINV, DET)
      DO I = 1, 3
        DO J = 1, 3
          FINV_T(I,J) = FINV(J,I)
        END DO
      END DO
      CALL MATMUL_3X3(FINV, SIG, VEC)
      CALL MATMUL_3X3(VEC, FINV_T, PK2)
      DO I = 1, 3
        DO J = 1, 3
          PK2(I,J) = AJ * PK2(I,J)
        END DO
      END DO
C
      RETURN
      END


C ============================================================
C  3×3 matrix utilities
C ============================================================
      SUBROUTINE INV3X3(A, AINV, DET)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3), AINV(3,3), DET
      DOUBLE PRECISION COF(3,3)
C
      COF(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
      COF(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
      COF(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
      COF(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
      COF(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
      COF(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
      COF(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
      COF(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
      COF(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C
      DET = A(1,1)*COF(1,1) + A(1,2)*COF(2,1) + A(1,3)*COF(3,1)
      IF(DET .EQ. 0.0D0) DET = 1.0D-20
C
      AINV(1,1) = COF(1,1)/DET
      AINV(2,1) = COF(2,1)/DET
      AINV(3,1) = COF(3,1)/DET
      AINV(1,2) = COF(1,2)/DET
      AINV(2,2) = COF(2,2)/DET
      AINV(3,2) = COF(3,2)/DET
      AINV(1,3) = COF(1,3)/DET
      AINV(2,3) = COF(2,3)/DET
      AINV(3,3) = COF(3,3)/DET
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
