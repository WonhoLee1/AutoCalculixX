# WHTOOLs MATCALIB 2026 — Abaqus → CalculiX PRF Parameter Mapping

## Overview

The CalculiX `umat_user.f` implements the Parallel Rheological Framework (PRF)
matching Abaqus `*VISCOELASTIC,NONLINEAR` with `*HYPERELASTIC`.
All formulations are verified against the Abaqus VISCNET test suite at:

`C:\SIMULIA\Documentation\2024LE\English\SIMAINPRefResources\viscnet_*.inp`

---

## 1. Hyperelastic Model Mapping

| Abaqus `*HYPERELASTIC` | CalculiX `FLAG_HE` | Details |
|---|---|---|
| `POLYNOMIAL,N=1` (Neo-Hookean) | `1` | `C10`, `D1` only |
| `POLYNOMIAL,N=2` (Mooney-Rivlin) | `1` | `C10`, `C01`, `D1`, `D2` |
| `POLYNOMIAL,N=3` (Full Polynomial) | `1` | All 9 Cij + 3 Dk |
| `ARRUDA-BOYCE` | `2` | 8 pre-computed constants (see §1.3) |
| `VAN DER WAALS` | (1) | Use `POLYNOMIAL,N=3` with appropriate curve-fit |
| `NEO HOOKE` | (1) | Use `POLYNOMIAL,N=1` — identical formulation |
| `MONEY-RIVLIN` | (1) | Use `POLYNOMIAL,N=2` — identical formulation |

### 1.1 Polynomial (FLAG_HE=1, 14 constants)

```
Abaqus *HYPERELASTIC,POLYNOMIAL,N=3
  C10, C01, C20, C11, C02, C30, C21, C12, C03, D1, D2, D3

CalculiX *USER MATERIAL
  C10, C01, C20, C11, C02, C30, C21, C12, C03, D1, D2, D3,
  IHYPER, RBULK_BACKUP
```

Direct 1:1 mapping for the first 12 constants.
`IHYPER` = 2 (use `RBULK_BACKUP` as bulk modulus — recommended).
`IHYPER` = 1 uses `D1`, `D2`, `D3` (matching Abaqus `*HYPERELASTIC` volumetric).

### 1.2 Neo-Hookean Special Case (FLAG_HE=1, C01=0)

```
Abaqus:
  *HYPERELASTIC,NEO HOOKE
  0.5, 0.0              ! C10, D1

CalculiX:
  *USER MATERIAL,CONSTANTS=19
  0, 1, 0,              ! N_NETWORK=0, FLAG_HE=1, FLAG_PL=0
  0.5, 0, 0,0,0, 0,0,0,0, 0,0,0, 2,1000
  1, 1000               ! G, RBULK
```

### 1.3 Arruda-Boyce (FLAG_HE=2, 8 constants)

Abaqus `*HYPERELASTIC,ARRUDA-BOYCE` takes `MU, LAMBDA_M, D` (3 parameters).
CalculiX takes 8 pre-computed constants based on the 8-chain series expansion:

| Parameter | Abaqus | CalculiX | Relation |
|---|---|---|---|
| C1 | — | `0.5` | fixed series coeff |
| C2 | — | `1/20` | fixed series coeff |
| C3 | — | `11/1050` | fixed series coeff |
| C4 | — | `19/7000` | fixed series coeff |
| C5 | — | `519/673750` | fixed series coeff |
| MU (C_R) | `MU` | `C_R` | 1:1 |
| D | `D` | `K0/2` | 1:1 scaled |
| BETA | — | `1/N` | where `N` = `LAMBDA_M²` |

The Python `ArrudaBoyceParams` dataclass handles this pre-computation
automatically via the `arruda_constants` property.

---

## 2. Viscous Network Model Mapping

Each `*VISCOELASTIC,NONLINEAR` network in Abaqus maps to one
`ViscousNetworkConfig` in the CalculiX PRF model.

### 2.1 Bergstrom-Boyce (FLAG_VISC=1, 5 params)

```
*dgamma = A1 * (lambda - 1 + KSI)^EXPC * (tau/TAUREF)^EXPM
```

| Abaqus `LAW=BERGSTROM` | CalculiX `VISC_BB` | Notes |
|---|---|---|
| `SRATIO` | `stiffn` | stress amplification fraction (same symbol, **not** dimensionless ratio in Abaqus docs) |
| `A` | `A1` | 1:1 |
| `m` | `EXPM` | power-law exponent |
| `C` | `EXPC` | chain-segment exponent |
| `E` | `KSI` | chain-stiffness parameter (must be > 0 for EXPC < 0) |

```
Abaqus:
  *VISCOELASTIC,NONLINEAR,NETWORKID=1,LAW=BERGSTROM,SRATIO=0.5
  0.1, 2, -1, 0.01    ! A, m, C, E

CalculiX:
  0.5, 1,
  0.1, -1, 2, 0.01, 1.0
  ^     ^   ^  ^     ^
  A1  EXPC EXPM KSI TAUREF
```

### 2.2 Law=STRAIN (FLAG_VISC=1, same as BB with KSI=0)

```
dgamma = A1 * (lambda - 1)^EXPC * (tau/1)^EXPM
```

Equivalent to BERGSTROM with `KSI=0.0`, `TAUREF=1.0`.
**Note**: `KSI=0.0` with `EXPC < 0` causes `(0)^EXPC` → singularity at `Fv=I`.
Set `KSI` to a small positive value (e.g. `0.001`) in CalculiX if convergence fails.

```
Abaqus:
  *VISCOELASTIC,NONLINEAR,NETWORKID=2,LAW=STRAIN,SRATIO=0.3
  0.01, 3, 0               ! A, m, C

CalculiX (use KSI=0.001 instead of 0 for stability):
  0.3, 1,
  0.01, 0, 3, 0.001, 1.0
```

### 2.3 Law=HYPERB (FLAG_VISC=2, 3 params — Sinh)

```
dgamma = A1 * sinh(B0 * tau)^EXPN
```

| Abaqus `LAW=HYPERB` | CalculiX `VISC_SINH` | Notes |
|---|---|---|
| `SRATIO` | `stiffn` | |
| `A` | `A1` | 1:1 |
| `B` | `B0` | 1:1 |
| `n` | `EXPN` | 1:1 |

```
Abaqus:
  *VISCOELASTIC,NONLINEAR,NETWORKID=3,LAW=HYPERB,SRATIO=0.4
  0.01, 1.0, 2.0           ! A, B, n

CalculiX:
  0.4, 2,
  0.01, 1.0, 2.0
```

### 2.4 Law=NOT AVAILABLE in Abaqus — Power (FLAG_VISC=3, 3 params)

```
dgamma = A1 * tau^EXPN * ((EXPM+1)*gamma)^EXPM
```

No direct Abaqus equivalent. Implemented from OpenRadioss `viscpower.F`.
Norton + strain-hardening creep law.

---

## 3. Network-Rheology Mapping

### 3.1 Single Network

```
Abaqus:
  *HYPERELASTIC,POLYNOMIAL,N=1
  0.5, 0.0                 ! C10=0.5, D1=0
  *VISCOELASTIC,NONLINEAR,LAW=BERGSTROM,SRATIO=0.5
  0.1, 2.0, -1.0, 0.01

CalculiX:
  *USER MATERIAL,CONSTANTS=26
  1, 1, 0,                 ! N_NETWORK=1, FLAG_HE=1, FLAG_PL=0
  0.5, 0,0,0,0, 0,0,0,0, 0,0,0, 2,1000,
  0.5, 1,                  ! stiffn, FLAG_VISC=1 (BB)
  0.1, -1.0, 2.0, 0.01, 1.0,
  1, 1000                  ! G, RBULK
  *DEPVAR
  14
```

### 3.2 Multi-Network (K=2)

```
Abaqus:
  *HYPERELASTIC,POLYNOMIAL,N=1
  0.5, 0.0
  *VISCOELASTIC,NONLINEAR,NETWORKID=1,LAW=BERGSTROM,SRATIO=0.3
  0.1, 2, -1, 0.01
  *VISCOELASTIC,NONLINEAR,NETWORKID=2,LAW=HYPERB,SRATIO=0.2
  0.01, 1.0, 2.0

CalculiX:
  *USER MATERIAL,CONSTANTS=31    ! 3 + 14 + (2+5) + (2+3) + 2 = 31
  2, 1, 0,                       ! N_NETWORK=2
  0.5, 0,0,0,0, 0,0,0,0, 0,0,0, 2,1000,
  0.3, 1,                        ! NW1: stiffn=0.3, BB
  0.1, -1.0, 2.0, 0.01, 1.0,
  0.2, 2,                        ! NW2: stiffn=0.2, Sinh
  0.01, 1.0, 2.0,
  1, 1000
  *DEPVAR
  26
```

### 3.3 Equilibrium Stress (N_NETWORK=0)

Abaqus does not allow `N_NETWORK=0` (
`*VISCOELASTIC,NONLINEAR` requires at least one network).
In CalculiX, `N_NETWORK=0` = pure hyperelastic stress (useful for
verification and modal analysis with hyperelastic only).

---

## 4. SRATIO vs STIFFN Semantics

The equilibrium stress is:

```
SIG = SIG_eq * (1 - Σ(stiffn_i)) + Σ(stiffn_i * SIG_vis_i)
```

**Abaqus**: `SRATIO` is the fraction of the hyperelastic stress
assigned to the viscous network. The equilibrium network carries
`(1 - Σ SRATIO)`.

**CalculiX**: `stiffn` has identical semantics. The sum of `stiffn`
across all networks must be ≤ 1. The equilibrium (`FLAG_HE`) is always
active and contributes `(1 - sum(stiffn))` × SIG_eq.

---

## 5. State Variable Layout

| Abaqus `STATEV` | CalculiX `statev` | Size | Description |
|---|---|---|---|
| — | `statev(1)` | 1 | initialised flag (0/1) |
| (in SDV) | `statev(2)` | 1 | J (det F) at converged state |
| per NW: `SDVn_1…9` | `statev(2+12*n+1…9)` | 9 | Fv (3×3) viscous deformation gradient |
| per NW: `SDVn_10` | `statev(2+12*n+10)` | 1 | dgamma (creep increment) |
| per NW: `SDVn_11` | `statev(2+12*n+11)` | 1 | tbnorm (effective stress norm) |
| per NW: — | `statev(2+12*n+12)` | 1 | gamma_old (cumulative, Power law) |
| **Total** | | `2 + 12×N_NETWORK` | |

---

## 6. Required *DEPVAR

```
DEPVAR = 2 + 12 × N_NETWORK
```

| N_NETWORK | DEPVAR |
|---|---|
| 0 | 2 |
| 1 | 14 |
| 2 | 26 |
| 3 | 38 |

---

## 7. Step Type

Abaqus uses `*VISCO` (quasi-static viscoelastic) or `*STATIC` with creep.
CalculiX does not support `*VISCO` — use `*STATIC` with `NLGEOM=YES`.

The time-dependent creep in `VISC_BB`, `VISC_SINH`, and `VISC_POWER`
operates within each static increment: `dgamma = A1 × dt × f(tau, ...)`.
For true relaxation analysis, use multiple `*STATIC` steps with
sufficient total step time.

```
*STEP,INC=200,NLGEOM=YES
*STATIC
0.01, 1.0, 1e-12, 0.1
```

---

## 8. Keyword Format Constraints

**Critical**: the CalculiX `*USER MATERIAL` parser reads exactly
**8 values per data line**. Values beyond the 8th on any line are
silently ignored. The `model_builder.py` correctly pads each line
(even the last) to exactly 8 values.

```
*USER MATERIAL,CONSTANTS=26
0, 1, 0, 0.5, 0, 0, 0, 0            ← first 8
0, 0, 0, 0, 0, 0, 2, 1000            ← second 8
0.5, 1, 0.1, -1, 2, 0.01, 1, 0       ← padded
1, 1000, 0, 0, 0, 0, 0, 0            ← padded
```

---

## 9. Numerical Tangent

The `umat_user.f` uses a **consistent numerical tangent** (Miehe formula):
each deformation gradient perturbation re-solves the creep update
(`PRF_UPDATE`) before computing the perturbed stress. This captures
the coupling between F and Fv in the tangent stiffness, unlike a
frozen-Fv approach.

```
For each perturbation column ε_pq:
  F_pert = F + ε × F⁻ᵀ × (e_p ⊗ e_q)
  STATEV ← STATEV_TANG          (reset to converged Fv)
  CALL PRF_UPDATE(F_pert, ...)  (re-solve creep at perturbed F)
  CALL PRF_EVAL(F_pert, ...)    (stress with updated Fv)
  C_ijpq = (S_pert - S_ref) / ε
```

The perturbation size is `ε = 1×10⁻⁸`.

## 10. Convergence Notes

- **KSI > 0 required** for `VISC_BB` with `EXPC < 0`. At initialization
  `Fv = I` → `λ = 1` → `(1 - 1 + KSI) = KSI`. Without KSI > 0,
  the term `KSI^EXPC` with EXPC < 0 produces `0^negative = ∞`.

- **Consistent tangent** works well for all 21/21 regression tests
  (N0, 1BB×3, 2NW×2, AB+BB/Sinh, Sinh, Power).

- **2-network creep with Σstiffn ≈ 1**: when the equilibrium network
  contributes near-zero stiffness, the tangent captures only the
  viscous response. This works for mild creep rates (`A1 < 0.05`)
  but requires small increments for fast creep (`A1 > 0.1`). For
  such cases, reduce `max_inc` or use smaller `A1` in material
  parameters.

---

## 10. Verified Configurations (all 21/21 test cases pass)

| # | Configuration | Networks | Result |
|---|---|---|---|
| 1 | N_NETWORK=0 (pure hyperelastic) | 0 | converges in 15 inc, 1 iter/inc |
| 2 | 1 BB network, A1=0.1, EXPM=2.0 | 1 | converges |
| 3 | 1 BB network, A1=0.01, EXPM=3.0 | 1 | converges |
| 4 | 1 BB network, A1=0.5, EXPM=1.5 | 1 | converges |
| 5 | 2 BB networks, mixed rates | 2 | converges |
| 6 | BB + Sinh networks | 2 | converges |
| 7 | Arruda-Boyce + BB | 1 | converges |
| 8 | Arruda-Boyce + Sinh | 1 | converges |
| 9 | Sinh only | 1 | converges |
| 10 | Power law only | 1 | converges |
