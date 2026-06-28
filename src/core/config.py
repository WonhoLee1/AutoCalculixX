# -*- coding: utf-8 -*-
# WHTOOLs MATCALIB 2026 — Material configuration module
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

BASE_DIR = Path(__file__).resolve().parent.parent.parent

# Use ccx_custom.exe which includes the hyperelastic + PRF UMAT
CALCULIX_EXE = r"D:\SOFTWARE\calculix_2.23_4win\ccx_custom.exe"
WORKSPACE_DIR = BASE_DIR / "workspace"


@dataclass
class ModalAnalysisConfig:
    job_name: str = "modal_job"
    num_modes: int = 10
    # Material
    E: float = 210000.0       # MPa
    nu: float = 0.3
    rho: float = 7.85e-9      # ton/mm³
    # Shell
    thickness: float = 2.0    # mm
    elset_name: str = "ALL_SHELLS"


@dataclass
class TrayGeometryConfig:
    length: float = 100.0
    width: float = 50.0
    mesh_size: float = 5.0
    elset_name: str = "TrayShell"


@dataclass
class HyperelasticConfig:
    """Polynomial hyperelastic material configuration.

    Corresponds to OpenRadioss LAW100 / CalculiX umat_user interface.

    Strain energy:
      W = sum Cij * (I1_bar-3)^i * (I2_bar-3)^j + sum Dk * (J-1)^{2k}

    Special cases:
      C10 only                     -> Neo-Hookean
      C10 + C01                    -> Mooney-Rivlin (2-param)
      All coefficients nonzero     -> Full Polynomial (N=3)
    """
    job_name: str = "hyperelastic_job"

    # Polynomial coefficients (deviatoric)
    C10: float = 0.5
    C01: float = 0.0
    C20: float = 0.0
    C11: float = 0.0
    C02: float = 0.0
    C30: float = 0.0
    C21: float = 0.0
    C12: float = 0.0
    C03: float = 0.0

    # Volumetric coefficients
    D1: float = 0.0
    D2: float = 0.0
    D3: float = 0.0

    # Volume mode: 1=D1/D2/D3 coefficients, 2=bulk modulus (recommended)
    IHYPER: int = 2
    RBULK: float = 1000.0    # bulk modulus (used when IHYPER=2)

    # Mesh / element
    elset_name: str = "ELSET1"
    element_type: str = "C3D8"

    # Step control (static NLGEOM)
    initial_inc: float = 0.05
    step_time: float = 1.0
    min_inc: float = 1e-8
    max_inc: float = 0.1

    @property
    def material_name(self) -> str:
        """Must start with 'USER' for umat_main dispatch."""
        return "USERHYPER"

    @property
    def constants(self) -> list:
        """Returns the 14 constants in elconloc order."""
        return [
            self.C10, self.C01, self.C20, self.C11, self.C02,
            self.C30, self.C21, self.C12, self.C03,
            self.D1,  self.D2,  self.D3,
            float(self.IHYPER), self.RBULK,
        ]


# ---------------------------------------------------------------------------
# PRF (Parallel Rheological Framework) — Multi-Network Foam UMAT
# ---------------------------------------------------------------------------

VISC_BB: int = 1    # Bergstrom-Boyce creep
VISC_SINH: int = 2  # Sinh-law creep
VISC_POWER: int = 3 # Power-law creep
VISC_PRONY: int = 4 # Prony series (linear viscoelastic)


@dataclass
class ViscousNetworkConfig:
    """Single viscous network within the PRF model.

    `stiffn` is the Abaqus SRATIO: fraction of total hyperelastic stress
    assigned to this network. Equilibrium carries (1 - sum(stiffn)).

    Creep equations matched to Abaqus VISCOELASTIC,NONLINEAR:

    BB (Bergstrom-Boyce, flag_visc=1): A1, EXPC, EXPM, KSI, TAUREF
      dgamma = A1 * (lambda - 1 + KSI)^EXPC * (tau/TAUREF)^EXPM
      Matches Abaqus LAW=BERGSTROM and LAW=STRAIN.
      **KSI must be > 0** (e.g. 0.001) to prevent singularity at Fv=I
      when EXPC < 0. At initialization Fv=I => lambda=1 => (1-1+KSI)=KSI.
      Abaqus default: KSI=0.01 for BERGSTROM, KSI=0.0 for STRAIN.

    Sinh (flag_visc=2): A1, B0, EXPN
      dgamma = A1 * sinh(B0 * tau)^EXPN
      Matches Abaqus LAW=HYPERB.

    Power (flag_visc=3): A1, EXPN, EXPM
      dgamma = A1 * ((EXPM+1)*gamma_old)^EXPM * tau^EXPN]^(1/(1+EXPM))
      Integrated Norton + strain-hardening. Matches OpenRadioss POWER.
      Does not match Abaqus LAW=STRAIN (use BB with EXPC=-n for that).

    Abaqus PRF parameter equivalents:
      *VISCOELASTIC,NONLINEAR,NETWORKID=1,SRATIO=R,LAW=STRAIN
        A, m, C  ->  VISC_BB: A1=A, EXPM=m, EXPC=C, KSI=0.0, TAUREF=1.0
      *VISCOELASTIC,NONLINEAR,NETWORKID=2,SRATIO=R,LAW=BERGSTROM
        A, m, C, E  ->  VISC_BB: A1=A, EXPM=m, EXPC=C, KSI=E, TAUREF=1.0
      *VISCOELASTIC,NONLINEAR,NETWORKID=3,SRATIO=R,LAW=HYPERB
        A, B, n  ->  VISC_SINH: A1=A, B0=B, EXPN=n
    """
    stiffn: float = 1.0        # SRATIO fraction (time-independent network)
    flag_visc: int = VISC_BB   # 1=BB, 2=Sinh, 3=Power

    # --- BB parameters ---
    A1: float = 0.001
    EXPC: float = -1.0
    EXPM: float = 2.0
    KSI: float = 0.001         # must be >0 for EXPC<0 (singularity at Fv=I)
    TAUREF: float = 1.0

    # --- Sinh parameters ---
    B0: float = 1.0

    # --- Power-law parameters ---
    EXPN: float = 2.0

    # --- Prony parameters (FLAG_VISC=4) ---
    G_PRONY: float = 0.5    # normalized shear relaxation modulus g_i
    K_PRONY: float = 0.0    # normalized bulk relaxation modulus k_i
    TAU_PRONY: float = 1.0  # relaxation time tau_i

    @property
    def nvisc(self) -> int:
        if self.flag_visc == VISC_BB:
            return 5
        elif self.flag_visc == VISC_SINH:
            return 3
        elif self.flag_visc == VISC_POWER:
            return 3
        elif self.flag_visc == VISC_PRONY:
            return 3    # g, k, tau
        return 0


@dataclass
class PolynomialHyperelasticParams:
    """Pr (deviatoric polynomial) + Qk (volumetric polynomial) params.
    Sent when PRFConfig.flag_he=1.
    """
    C10: float = 0.5
    C01: float = 0.0
    C20: float = 0.0
    C11: float = 0.0
    C02: float = 0.0
    C30: float = 0.0
    C21: float = 0.0
    C12: float = 0.0
    C03: float = 0.0

    D1: float = 0.0
    D2: float = 0.0
    D3: float = 0.0

    IHYPER: int = 2         # 1=D-coeffs, 2=RBULK
    RBULK_HE: float = 1000.0


@dataclass
class ArrudaBoyceParams:
    """Arruda-Boyce 8-chain parameters.
    Sent when PRFConfig.flag_he=2.

    The 8 constants passed to *USER MATERIAL are pre-computed from
    simplified inputs (N, C_R, K0):
      [C1=0.5, C2=0.05, C3=11/1050, C4=19/7000, C5=519/673750,
       MU=C_R, D=K0/2, BETA=1/N]
    Note: N here is locking stretch LAMBDA_M (matching Abaqus convention).
    The actual number of chain segments is N² = LAMBDA_M².
    """
    N: float = 2.5          # locking stretch LAMBDA_M (number of chain segments = N²)
    C_R: float = 0.5        # rubbery modulus MU
    K0: float = 1000.0      # bulk modulus

    @property
    def arruda_constants(self) -> list:
        C1 = 0.5
        C2 = 1.0 / 20.0
        C3 = 11.0 / 1050.0
        C4 = 19.0 / 7000.0
        C5 = 519.0 / 673750.0
        MU = self.C_R
        D  = self.K0 / 2.0
        BETA = 1.0 / (self.N * self.N)
        return [C1, C2, C3, C4, C5, MU, D, BETA]


@dataclass
class YeohParams:
    """Yeoh (Reduced Polynomial 3rd order) hyperelastic params.
    Sent when PRFConfig.flag_he=3.

    W_dev = C10*(I1b-3) + C20*(I1b-3)^2 + C30*(I1b-3)^3
    U(J)  = D1*(J-1)^2 + D2*(J-1)^4 + D3*(J-1)^6
    """
    C10: float = 0.5
    C20: float = 0.0
    C30: float = 0.0
    D1: float = 0.0
    D2: float = 0.0
    D3: float = 0.0

    @property
    def constants(self) -> list:
        return [self.C10, self.C20, self.C30, self.D1, self.D2, self.D3]


@dataclass
class GentParams:
    """Gent hyperelastic (locking model) params.
    Sent when PRFConfig.flag_he=4.

    W_dev = -(MU/2) * Jm * ln(1 - (I1b-3)/Jm)
    U(J)  = D*(J-1)^2   (via dU/dJ = D*(J - 1/J))

    MU  = initial shear modulus
    Jm  = locking stretch (max I1b-3 before stiffening)
    D   = volumetric parameter (same as Arruda-Boyce D)
    """
    MU: float = 1.0
    Jm: float = 10.0
    D: float = 500.0

    @property
    def constants(self) -> list:
        return [self.MU, self.Jm, self.D]


@dataclass
class OgdenParams:
    """Ogden (N=3) hyperelastic params.
    Sent when PRFConfig.flag_he=5.

    W_dev = sum_{k=1..3} (2*mu_k/alpha_k^2) *
            (lam1b^alpha_k + lam2b^alpha_k + lam3b^alpha_k - 3)
    U(J)  = D1*(J-1)^2 + D2*(J-1)^4 + D3*(J-1)^6

    When alpha_k -> 0, term degenerates to mu_k * ln(lambda).
    """
    mu1: float = 1.0
    alpha1: float = 2.0
    mu2: float = 0.0
    alpha2: float = 0.0
    mu3: float = 0.0
    alpha3: float = 0.0
    D1: float = 0.0
    D2: float = 0.0
    D3: float = 0.0

    @property
    def constants(self) -> list:
        return [self.mu1, self.alpha1, self.mu2, self.alpha2,
                self.mu3, self.alpha3, self.D1, self.D2, self.D3]


@dataclass
class PRFConfig:
    """PRF (Parallel Rheological Framework) material configuration.

    Material constant layout sent to *USER MATERIAL:
      [N_NETWORK, FLAG_HE, FLAG_PL, <HE_params>, <NW1>, ..., <NWn>, G, RBULK]
    """
    job_name: str = "prf_job"

    # Hyperelastic model type
    flag_he: int = 1                       # 1=polynomial, 2=Arruda-Boyce,
                                           # 3=Yeoh, 4=Gent, 5=Ogden
    he_poly: PolynomialHyperelasticParams = field(
        default_factory=PolynomialHyperelasticParams)
    he_ab: ArrudaBoyceParams = field(
        default_factory=ArrudaBoyceParams)
    he_yeoh: YeohParams = field(
        default_factory=YeohParams)
    he_gent: GentParams = field(
        default_factory=GentParams)
    he_ogden: OgdenParams = field(
        default_factory=OgdenParams)

    # Viscous networks (list of ViscousNetworkConfig)
    networks: List[ViscousNetworkConfig] = field(default_factory=lambda: [
        ViscousNetworkConfig(),
    ])

    # Wave-speed check values (not used in stress computation)
    G: float = 1.0
    RBULK_FINAL: float = 1000.0

    # Mesh / element
    elset_name: str = "ELSET1"
    element_type: str = "C3D8"

    # Step control (static NLGEOM with creep)
    initial_inc: float = 0.01
    step_time: float = 1.0
    min_inc: float = 1e-12
    max_inc: float = 0.1
    max_disp: float = 0.3       # max displacement for multi-step plots

    @property
    def n_networks(self) -> int:
        return len(self.networks)

    @property
    def material_name(self) -> str:
        return "USERPRF"

    @property
    def nconstants(self) -> int:
        """Total number of *USER MATERIAL constants."""
        n = 3  # N_NET, FLAG_HE, FLAG_PL
        if self.flag_he == 1:
            n += 14  # polynomial
        elif self.flag_he == 2:
            n += 8   # Arruda-Boyce
        elif self.flag_he == 3:
            n += 6   # Yeoh
        elif self.flag_he == 4:
            n += 3   # Gent
        elif self.flag_he == 5:
            n += 9   # Ogden
        n += sum(2 + nw.nvisc for nw in self.networks)  # network blocks
        n += 2  # G, RBULK
        return n

    @property
    def depvar(self) -> int:
        """Number of state variables = 2 + 12 * N_NETWORK."""
        return 2 + 12 * self.n_networks

    @property
    def constants(self) -> list:
        c: list = [float(self.n_networks), float(self.flag_he), 0.0]  # FLAG_PL=0

        if self.flag_he == 1:
            p = self.he_poly
            c += [p.C10, p.C01, p.C20, p.C11, p.C02,
                  p.C30, p.C21, p.C12, p.C03,
                  p.D1,  p.D2,  p.D3,
                  float(p.IHYPER), p.RBULK_HE]
        elif self.flag_he == 2:
            a = self.he_ab
            c += a.arruda_constants
        elif self.flag_he == 3:
            c += self.he_yeoh.constants
        elif self.flag_he == 4:
            c += self.he_gent.constants
        elif self.flag_he == 5:
            c += self.he_ogden.constants

        for nw in self.networks:
            c += [nw.stiffn, float(nw.flag_visc)]
            if nw.flag_visc == VISC_BB:
                c += [nw.A1, nw.EXPC, nw.EXPM, nw.KSI, nw.TAUREF]
            elif nw.flag_visc == VISC_SINH:
                c += [nw.A1, nw.B0, nw.EXPN]
            elif nw.flag_visc == VISC_POWER:
                c += [nw.A1, nw.EXPN, nw.EXPM]
            elif nw.flag_visc == VISC_PRONY:
                c += [nw.G_PRONY, nw.K_PRONY, nw.TAU_PRONY]

        c += [self.G, self.RBULK_FINAL]
        return c
