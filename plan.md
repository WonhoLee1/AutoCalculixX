# AutoCalculix - LAW100 Hyperelastic UMAT 설계

## 0. 개요

**목표**: OpenRadioss의 `/MAT/LAW100` (Polynomial / Neo-Hookean / Mooney-Rivlin hyperelastic 재료 모델)을 CalculiX/Abaqus UMAT 사용자 서브루틴으로 포팅하여 기존 AutoCalculix 파이프라인에 통합.

**범위**: 
- Polynomial strain energy (범용) + Neo-Hookean (특수 케이스)
- 단일 요소 검증 + 기존 모달 파이프라인 연계
- Fortran ISO_C_BINDING 없이 순수 Fortran 77/90 UMAT 인터페이스

---

## 1. UMAT vs LAW100 인터페이스 차이 분석

| 항목 | OpenRadioss LAW100 | CalculiX UMAT |
|---|---|---|
| **입력** | Deformation gradient `F` (적분형) | Strain increment `DSTRAN` (증분형) |
| **출력** | Cauchy stress `σ` (full tensor) | Cauchy(=Kirchhoff in finite) `STRESS` (Voigt 6) |
| **Tangent** | 없음 (explicit dynamics) | `DDSDDE` (6×6 Jacobian) — implicit |
| **상태 변수** | `UVAR(NEL, NUVAR)` (array) | `STATEV(NSTATV)` (1D vector) |
| **재료 상수** | `UPARAM(NUPARAM)` (1D) | `PROPS(NPROPS)` (1D) |
| **시간** | `TIME`, `TIMESTEP` | `DTIME` only |
| **요소 루프** | Solver 외부 루프 (NEL batch) | UMAT 내부 단일 적분점 x 1회 호출 |
| **차원** | 3D solids only | 3D / plane strain / axisymmetric |
| **방향** | Global frame (corotational 따로) | Global frame (자체 처리) |
| **열 팽창** | `fth = I + αΔT` (직접 분리) | `DSTRAN`에서 solver가 분리 (UMAT은 기계적 변형만) |

**핵심 변환**: OpenRadioss의 `F` 기반 적분형 → CalculiX의 증분형 `DSTRAN` 기반
- 방법: UMAT 진입 시 `DFGRD0`(t=0 기준 F)와 `DFGRD1`(현재 F)를 사용 → **적분형으로 직접 계산 가능**
- CalculiX/Abaqus UMAT은 `DFGRD0`, `DFGRD1` 배열을 통해 F를 직접 받을 수 있음

---

## 2. 수학적 모델

### 2.1 Polynomial Strain Energy (범용)

```
W = Σ_{i+j=1}^{3} Cij (Ī₁ - 3)^i (Ī₂ - 3)^j
  + Σ_{k=1}^{3} Dk (J - 1)^{2k}

Ī₁ = J^{-2/3} · I₁        (deviatoric 1st invariant)
Ī₂ = J^{-4/3} · I₂        (deviatoric 2nd invariant)
J  = det(F)                (volume ratio)
```

### 2.2 특수 케이스

| 모델 | 조건 | W |
|---|---|---|
| **Neo-Hookean** | C10 ≠ 0, 나머지 0 | `W = C10(Ī₁ - 3) + D₁(J - 1)²` |
| **Mooney-Rivlin** | C10, C01 ≠ 0, 나머지 0 | `W = C10(Ī₁ - 3) + C01(Ī₂ - 3) + D₁(J - 1)²` |
| **Full Polynomial** | 모든 Cij | 위 일반식 |

### 2.3 Cauchy Stress (Kirchhoff → Cauchy)

```
σ = (2/J) [ (∂W/∂Ī₁ + ∂W/∂Ī₂ · Ī₁) · b̄ − ∂W/∂Ī₂ · b̄² ] 
  − (2/3J)(Ī₁·∂W/∂Ī₁ + 2Ī₂·∂W/∂Ī₂)·I + ∂W/∂J · I
```

여기서 `b̄ = J^{-2/3} · F·F^T` (deviatoric left Cauchy-Green tensor)

### 2.4 Tangent Modulus (DDSDDE) — Implicit Essential

6×6 consistent tangent를 analytical하게 계산 (또는 수치 미분 fallback):

```
ℂ = ℂ_vol + ℂ_dev

ℂ_vol_ijkl = (J·∂²W/∂J² + ∂W/∂J) δ_ij δ_kl - 2·∂W/∂J · 𝕀_ijkl
ℂ_dev     = (4/3)∂W/∂Ī₁ · ... (b̄ deviatoric projection)
```

→ 1차 구현에서는 **수치적 tangent (perturbation)** 사용 후 analytical로 교체

---

## 3. 파일 구조

```
D:\PythonCodeStudy\AutoCalculix\
├── src/
│   ├── umat/                          # NEW - Fortran UMAT 소스
│   │   ├── umat_hyperelastic.f        # 메인 UMAT (CalculiX 인터페이스)
│   │   ├── polynomial_stress.f        # POLYSTRESS2 포팅 (Polynomial 응력)
│   │   ├── polynomial_tangent.f       # Analytical tangent (2단계)
│   │   └── make_umat.bat             # Windows gfortran 빌드 스크립트
│   │
│   ├── core/
│   │   ├── config.py                  # 수정: HyperelasticConfig dataclass 추가
│   │   ├── model_builder.py           # 수정: *USER MATERIAL, *DEPVAR INP 생성
│   │   └── solver.py                  # 수정: UMAT DLL 로드 경로 전달
│   │
│   ├── tests/                         # NEW - 검증
│   │   ├── test_single_element.py     # 단일 요소: uniaxial / shear / volumetric
│   │   └── fixtures/                  # 참조 결과 (Abaqus 검증 등)
│   │
│   └── pipeline.py                    # 수정: hyperelastic 모드 추가
│
├── workspace/                         # INP + 결과 (gitignore)
└── plan.md                            # 이 파일
```

### 3.1 Fortran 소스 상세

#### `umat_hyperelastic.f` — 메인 UMAT

```fortran
SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,
     1 RPL, DDSDDT, DRPLDE, DRPLDT,
     2 STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME,
     3 NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,
     4 CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, JSTEP, KINC)
! ...
! PROPS 맵핑:
!   PROPS(1)  = IHYPER  (1=Polynomial)
!   PROPS(2)  = C10
!   PROPS(3)  = C01
!   PROPS(4:10) = C20, C11, C02, C30, C21, C12, C03
!   PROPS(11) = D1
!   PROPS(12) = D2
!   PROPS(13) = D3
!   PROPS(14) = K (bulk modulus fallback, D1=0 시 사용)
```

OpenRadioss의 `sigeps100.F90`에서 핵심 로직만 추출:
1. `DFGRD1`로 `F` 구성 (3×3 → 2D는 패딩)
2. `POLYSTRESS` 호출 → Cauchy stress
3. Voigt 변환 (3×3 tensor → 6-component vector)
4. Tangent 계산 (수치 미분 1차)
5. `STATEV` 업데이트 (필요시)

#### `polynomial_stress.f` — 응력 계산 (순수 함수)

- [sigpoly.F](https://github.com/OpenRadioss/OpenRadioss/blob/main/engine/source/materials/mat/mat100/sigpoly.F)의 `POLYSTRESS` 서브루틴을 거의 그대로 포팅
- `MATB`(left Cauchy-Green) → `SIG`(Cauchy stress)
- 종속성 제거: `implicit_f.inc` → 명시적 선언, `EM20` → `1.0D-20`

### 3.2 Python 측 변경

#### `config.py` 추가

```python
@dataclass
class HyperelasticConfig:
    """Polynomial hyperelastic material model configuration."""
    model_type: str = "polynomial"  # "polynomial" | "neo_hookean" | "mooney_rivlin"
    C10: float = 0.0
    C01: float = 0.0
    C20: float = 0.0
    C11: float = 0.0
    C02: float = 0.0
    C30: float = 0.0
    C21: float = 0.0
    C12: float = 0.0
    C03: float = 0.0
    D1: float = 0.0   # incompressibility
    D2: float = 0.0
    D3: float = 0.0
    rho: float = 1.0e-9
```

**기본값 예시**:
- Neo-Hookean: `C10 = μ/2 = G/2` (e.g., C10=0.5, D1=0.01)
- Mooney-Rivlin: `C10, C01` 둘 다 (e.g., C10=0.3, C01=0.1, D1=0.01)

#### `model_builder.py` 확장

```python
def _write_user_material(self, f, config: HyperelasticConfig):
    """*USER MATERIAL + *DEPVAR + *DENSITY INP 블록 생성"""
    nprops = 14
    nstatv = 1  # 최소 상태 변수
    f.write(f"*USER MATERIAL, CONSTANTS={nprops}\n")
    vals = [config.model_type_code, config.C10, config.C01, ...]
    f.write(", ".join(str(v) for v in vals) + "\n")
    f.write(f"*DEPVAR\n{nstatv}\n")
    f.write(f"*DENSITY\n{config.rho}\n")
```

---

## 4. 검증 전략

### 4.1 단일 요소 검증 (C3D8 / C3D20R)

| 테스트 | 변형 모드 | 검증 방법 |
|---|---|---|
| Uniaxial tension | ε₁₁ = 0.5 | Analytical σ = 2C10(λ - λ⁻²) + 2C01(1 - λ⁻³) (Mooney-Rivlin) |
| Equibiaxial tension | ε₁₁ = ε₂₂ = 0.3 | Analytical formula 비교 |
| Simple shear | γ₁₂ = 0.5 | µ = 2(C10 + C01) → τ ≈ µγ (small strain) |
| Volumetric | J = 1.05 | σ_hydro = 2D₁(J-1) (Neo-Hookean) |

### 4.2 검증 프로세스

1. Python으로 단일 요소 INP 생성 (`*USER MATERIAL` 포함)
2. CalculiX 실행 → UMAT 로드 → 결과 확인
3. `test_single_element.py`로 자동화

---

## 5. 구현 순서

### Phase 1: Core Fortran UMAT (순수 함수)
1. `polynomial_stress.f` — POLYSTRESS 포팅 (상태 없음, 입력→출력)
2. `umat_hyperelastic.f` — UMAT 인터페이스 + 수치 tangent
3. 빌드 스크립트 (`make_umat.bat` — gfortran으로 DLL 생성)

### Phase 2: Python 통합
4. `config.py` — `HyperelasticConfig` dataclass
5. `model_builder.py` — `*USER MATERIAL` INP 생성 확장
6. `solver.py` — UMAT 환경변수/경로 설정

### Phase 3: 검증
7. `test_single_element.py` — 단일 요소 4종 테스트
8. 실행 및 결과 리포트

### Phase 4 (Optional): 개선
9. Analytical tangent (DDSDDE)
10. Mooney-Rivlin / Arruda-Boyce 확장
11. Modal analysis 통합 검증

---

## 6. 상세 기술 결정

### 6.1 Deformation gradient 처리

CalculiX UMAT은 `DFGRD0`(step 시작 F)와 `DFGRD1`(현재 F)를 제공하므로 **적분형 접근법** 가능:

```fortran
! F = DFGRD1
! b = F·F^T (left Cauchy-Green tensor)
! J = det(F)
```

이 접근법은 OpenRadioss의 적분형과 동일하여 검증이 용이함.

### 6.2 Voigt 변환 규칙

Abaqus/CalculiX Voigt order:
- Direct: 11, 22, 33, 12, 13, 23
- Shear strain: **engineering shear** γ (DSTRAN은 engineering)

OpenRadioss → CalculiX:
```fortran
STRESS(1) = SIG(1,1)  ! σ₁₁
STRESS(2) = SIG(2,2)  ! σ₂₂
STRESS(3) = SIG(3,3)  ! σ₃₃
STRESS(4) = SIG(1,2)  ! σ₁₂
STRESS(5) = SIG(1,3)  ! σ₁₃ (mirror)
STRESS(6) = SIG(2,3)  ! σ₂₃
```

**주의**: OpenRadioss는 σ₂₃, σ₃₁ 순서 → CalculiX는 σ₁₃, σ₂₃ 순서

### 6.3 Tangent 수치 미분 (1차 구현)

Perturbation approach:
```fortran
eps = 1.0D-6 * (1.0D0 + norm(STRAN))
DO i = 1, NTENS
    DSTRAN_pert = DSTRAN
    DSTRAN_pert(i) = DSTRAN_pert(i) + eps
    ! 재계산 → stress_pert
    DDSDDE(:,i) = (stress_pert - stress) / eps
END DO
```

비용: UMAT 호출 1회당 추가로 6회의 stress evaluation (총 ~7회). 
단일 요소 검증 단계에서는 충분.

### 6.4 빌드 방법 (Windows)

```batch
:: make_umat.bat
gfortran -c -O2 -fPIC polynomial_stress.f
gfortran -c -O2 -fPIC umat_hyperelastic.f
gfortran -shared -o umat_hyperelastic.dll polynomial_stress.o umat_hyperelastic.o
```

CalculiX는 환경변수 `CCX_UMAT_PATH` 또는 INP의 `*USER MATERIAL`에 DLL 경로를 지정.

---

## 7. 위험 / 제한사항

| 위험 | 영향 | 완화 |
|---|---|---|
| Windows gfortran ABI 이슈 | CalculiX가 다른 컴파일러로 빌드됨 | CalculiX 자체 빌드 시 동일 gfortran 사용 |
| 2D 요소 (shell) 미지원 | UMAT은 3D solids만 | 필요시 UMAT 외 shell section으로 linear elastic fallback |
| 재료 불안정성 (Drucker) | DDSDDE 비양정 → 수렴 실패 | D1 파라미터로 충분한 체적 강성 확보 |
| OpenRadioss 종속성 | `constant_mod`, `precision_mod` | 상수 직접 정의 (EM20=1e-20, ONE=1.0D0 등) |
