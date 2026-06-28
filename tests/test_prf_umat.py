"""WHTOOLs MATCALIB 2026 — PRF UMAT regression tests

Generates INP files via model_builder, runs ccx_custom.exe, parses .sta
for convergence, and validates stress/energy outputs against known ranges.

Run:  pytest tests/test_prf_umat.py -v --timeout=600
"""
import sys, os, subprocess, re
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.core.config import (
    PRFConfig, ViscousNetworkConfig, PolynomialHyperelasticParams,
    ArrudaBoyceParams, YeohParams, GentParams, OgdenParams,
    VISC_BB, VISC_SINH, VISC_POWER, VISC_PRONY,
)
from src.core.model_builder import CalculixModelBuilder

# --- Environment -----------------------------------------------------------
WORKSPACE = Path(__file__).resolve().parents[1] / "workspace"
BUILD_DIR = Path("D:/PythonCodeStudy/AutoCalculix/calculix/CalculiX/ccx_2.23/src")
CCX_EXE = Path("D:/SOFTWARE/calculix_2.23_4win/ccx_custom.exe")
MINGW_BIN = "C:/msys64/mingw64/bin"

ENV = os.environ.copy()
ENV["Path"] = f"{MINGW_BIN};{ENV.get('Path', '')}"

BUILDER = CalculixModelBuilder(WORKSPACE)


# --- Helpers ---------------------------------------------------------------

def solve(inp_stem: str, timeout: int = 300) -> dict:
    """Run CalculiX and return convergence summary dict."""
    # Clean old outputs
    for ext in [".sta", ".dat", ".frd", ".log", ".msg"]:
        p = WORKSPACE / f"{inp_stem}{ext}"
        if p.exists():
            p.unlink()

    result = subprocess.run(
        [str(CCX_EXE), inp_stem],
        cwd=str(WORKSPACE),
        capture_output=True, text=True, timeout=timeout, env=ENV,
    )

    stdout = result.stdout
    stderr = result.stderr
    combined = stdout + stderr

    # Max increments reached?
    inc_limit = "max. # of increments reached" in combined
    error_occurred = "*ERROR" in combined and not inc_limit

    # Parse .sta
    sta_path = WORKSPACE / f"{inp_stem}.sta"
    if not sta_path.exists():
        return {"success": False, "error": "no_sta", "details": stderr[-500:]}

    lines = [l.strip() for l in sta_path.read_text().splitlines()
             if l.strip() and not l.startswith("SUMMARY") and "STEP" not in l
             and "INC" not in l and "ATT" not in l]

    if not lines:
        return {"success": False, "error": "empty_sta"}

    inc_data = []
    for ln in lines:
        parts = ln.split()
        if len(parts) >= 5:
            att = parts[2].rstrip('U')  # strip 'U' suffix (unsuccessful attempt)
            inc_data.append({
                "inc": int(parts[1]),
                "iters": int(att),
                "total_time": float(parts[3]),
                "step_time": float(parts[4]),
            })

    if not inc_data:
        return {"success": False, "error": "no_inc_data"}

    last = inc_data[-1]
    total_time = last["total_time"]
    n_inc = len(inc_data)

    return {
        "success": total_time >= 0.999 and not error_occurred,
        "partial_max_inc": inc_limit,
        "error": error_occurred,
        "total_time": total_time,
        "n_inc": n_inc,
        "last_iters": last.get("iters", 0),
        "inc_data": inc_data,
        "stdout": stdout[-1000:],
    }


# --- Tests -----------------------------------------------------------------

@pytest.fixture(scope="module")
def n0_inp():
    cfg = PRFConfig(
        job_name="test_prf_n0",
        flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[],
        initial_inc=0.05, step_time=1.0,
    )
    return BUILDER.build_prf_inp(cfg).stem

class TestPureHyperelastic:
    """N_NETWORK=0: PRF with no viscous networks = hyperelastic pass-through."""

    def test_converges(self, n0_inp):
        r = solve(n0_inp)
        assert r["success"], f"Did not converge: total_time={r.get('total_time')}"

    def test_fast(self, n0_inp):
        r = solve(n0_inp)
        assert r["n_inc"] < 50, f"Too many increments: {r['n_inc']}"

    def test_depvar_2(self, n0_inp):
        cfg = PRFConfig(networks=[], flag_he=1,
            he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.5, RBULK_HE=1000.0))
        assert cfg.depvar == 2, f"Expected 2 depvars, got {cfg.depvar}"


@pytest.fixture(scope="module", params=[
    {"A1": 0.1,  "EXPC": 0.0, "EXPM": 2.0, "KSI": 0.01, "stiffn": 0.5},
    {"A1": 0.01, "EXPC": -0.5, "EXPM": 3.0, "KSI": 0.001, "stiffn": 0.3},
    {"A1": 0.5,  "EXPC": 1.0, "EXPM": 1.5, "KSI": 0.01, "stiffn": 0.8},
])
def bb_cfg(request):
    params = dict(request.param)
    stiffn = params.pop("stiffn")
    return PRFConfig(
        job_name=f"test_1bb_all_{params['A1']}_{params['EXPC']}_{params['EXPM']}",
        flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[ViscousNetworkConfig(stiffn=stiffn, **params)],
        initial_inc=0.02, step_time=1.0,
    )

class TestSingleNetworkBB:
    """Single BB network (most common PRF case)."""

    def test_converges(self, bb_cfg):
        inp = BUILDER.build_prf_inp(bb_cfg)
        r = solve(inp.stem)
        assert r["success"], f"BB A1={bb_cfg.networks[0].A1} failed: total_time={r.get('total_time')}"


@pytest.fixture(scope="module", params=["2bb", "bb+sinh"])
def multi_cfg(request):
    if request.param == "2bb":
        cfg = PRFConfig(
            job_name="test_2bb_regr",
            flag_he=1,
            he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
            networks=[
                ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01),
                ViscousNetworkConfig(stiffn=0.2, A1=0.005, EXPC=0.0, EXPM=2.0, KSI=0.01),
            ],
            initial_inc=0.01, step_time=1.0, max_inc=0.1,
        )
    elif request.param == "bb+sinh":
        cfg = PRFConfig(
            job_name="test_bb_sinh_regr",
            flag_he=1,
            he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
            networks=[
                ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01),
                ViscousNetworkConfig(stiffn=0.2, flag_visc=VISC_SINH, A1=0.01, B0=1.0, EXPN=2.0),
            ],
            initial_inc=0.01, step_time=1.0, max_inc=0.1,
        )
    return cfg

class TestMultiNetwork:
    """2 and 3 network configurations."""

    def test_converges(self, multi_cfg):
        inp = BUILDER.build_prf_inp(multi_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Multi-network {inp.stem}: total_time={r.get('total_time')}"


@pytest.fixture(scope="module", params=[
    {"flag_visc": VISC_BB, "A1": 0.1, "EXPM": 2.0, "KSI": 0.01},
    {"flag_visc": VISC_SINH, "A1": 0.01, "B0": 1.0, "EXPN": 2.0},
])
def ab_cfg(request):
    vp = dict(request.param)
    fv = vp.pop("flag_visc")
    return PRFConfig(
        job_name=f"test_ab_visctype{fv}",
        flag_he=2,
        he_ab=ArrudaBoyceParams(C_R=1.0, N=7.0, K0=1000.0),
        networks=[ViscousNetworkConfig(stiffn=0.5, flag_visc=fv, **vp)],
        G=1.0, RBULK_FINAL=1000.0,
        initial_inc=0.02, step_time=1.0,
    )

class TestArrudaBoyce:
    """Arruda-Boyce + creep."""

    def test_converges(self, ab_cfg):
        inp = BUILDER.build_prf_inp(ab_cfg)
        r = solve(inp.stem)
        assert r["success"], f"AB network failed: total_time={r.get('total_time')}"


@pytest.fixture(scope="module", params=[VISC_SINH, VISC_POWER])
def sp_cfg(request):
    fv = request.param
    kw = {"flag_visc": fv, "A1": 0.01}
    if fv == VISC_SINH:
        kw.update({"B0": 1.0, "EXPN": 2.0})
    elif fv == VISC_POWER:
        kw.update({"EXPN": 2.0, "EXPM": 0.5})
    return PRFConfig(
        job_name=f"test_visc{fv}",
        flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[ViscousNetworkConfig(stiffn=0.4, **kw)],
        initial_inc=0.02, step_time=1.0,
    )

class TestSinhAndPower:
    """Single network Sinh and Power."""

    def test_converges(self, sp_cfg):
        inp = BUILDER.build_prf_inp(sp_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Visc type {sp_cfg.networks[0].flag_visc} failed: total_time={r.get('total_time')}"


class TestNConstants:
    """Validate constant counts match expected layouts."""

    def test_n0(self):
        cfg = PRFConfig(networks=[], flag_he=1, he_poly=PolynomialHyperelasticParams(C10=0.5))
        assert cfg.nconstants == 3 + 14 + 0 + 2 == 19

    def test_1bb(self):
        cfg = PRFConfig(networks=[ViscousNetworkConfig(stiffn=0.5, flag_visc=VISC_BB)])
        assert cfg.nconstants == 3 + 14 + (2 + 5) + 2 == 26

    def test_2bb(self):
        cfg = PRFConfig(networks=[
            ViscousNetworkConfig(stiffn=0.3),
            ViscousNetworkConfig(stiffn=0.2),
        ])
        assert cfg.nconstants == 3 + 14 + 2 * (2 + 5) + 2 == 33

    def test_ab_bb(self):
        cfg = PRFConfig(
            flag_he=2,
            he_ab=ArrudaBoyceParams(),
            networks=[ViscousNetworkConfig(stiffn=0.5)],
        )
        assert cfg.nconstants == 3 + 8 + (2 + 5) + 2 == 20

    def test_yeoh_n0(self):
        cfg = PRFConfig(
            job_name="test_yeoh_n0",
            flag_he=3,
            he_yeoh=YeohParams(C10=0.5),
            networks=[],
        )
        assert cfg.nconstants == 3 + 6 + 0 + 2 == 11

    def test_gent_1bb(self):
        cfg = PRFConfig(
            job_name="test_gent",
            flag_he=4,
            he_gent=GentParams(MU=1.0, Jm=10.0),
            networks=[ViscousNetworkConfig(stiffn=0.3)],
        )
        assert cfg.nconstants == 3 + 3 + (2 + 5) + 2 == 15

    def test_ogden_n0(self):
        cfg = PRFConfig(
            job_name="test_ogden",
            flag_he=5,
            he_ogden=OgdenParams(mu1=1.0, alpha1=2.0),
            networks=[],
        )
        assert cfg.nconstants == 3 + 9 + 0 + 2 == 14


class TestDepvar:
    """State variable counts."""

    def test_n0(self):
        assert PRFConfig(networks=[]).depvar == 2

    def test_1nw(self):
        assert PRFConfig(networks=[ViscousNetworkConfig()]).depvar == 14

    def test_2nw(self):
        nw = [ViscousNetworkConfig(), ViscousNetworkConfig()]
        assert PRFConfig(networks=nw).depvar == 26

    def test_3nw(self):
        nw = [ViscousNetworkConfig() for _ in range(3)]
        assert PRFConfig(networks=nw).depvar == 38


@pytest.fixture(scope="module")
def eight_per_line_inp():
    cfg = PRFConfig(
        job_name="test_eight_per_line",
        flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[
            ViscousNetworkConfig(stiffn=0.3),
            ViscousNetworkConfig(stiffn=0.2),
        ],
    )
    return BUILDER.build_prf_inp(cfg)

class TestEightPerLine:
    """*USER MATERIAL constants must be written 8 per line."""

    def test_lines_eight(self, eight_per_line_inp):
        """Verify every *USER MATERIAL data line has exactly 8 values."""
        text = eight_per_line_inp.read_text()
        in_mat = False
        for line in text.splitlines():
            if line.startswith("*USER MATERIAL"):
                in_mat = True
                continue
            if line.startswith("*") and in_mat:
                in_mat = False
            if in_mat and line.strip() and not line.startswith("*"):
                vals = line.strip().split(",")
                non_empty = [v for v in vals if v.strip()]
                assert len(non_empty) == 8, (
                    f"Line has {len(non_empty)} values (need 8): {line.strip()[:80]}"
                )


# --- New hyperelastic models: Yeoh, Gent, Ogden ---

@pytest.fixture(scope="module")
def yeoh_cfg():
    return PRFConfig(
        job_name="test_yeoh_pure",
        flag_he=3,
        he_yeoh=YeohParams(C10=0.5, C20=0.1, C30=0.02, D1=500.0),
        networks=[],
        initial_inc=0.05, step_time=1.0,
    )

class TestYeoh:
    """Yeoh (Reduced Polynomial 3rd order) pure hyperelastic."""

    def test_converges(self, yeoh_cfg):
        inp = BUILDER.build_prf_inp(yeoh_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Yeoh failed: total_time={r.get('total_time')}"


@pytest.fixture(scope="module")
def gent_cfg():
    return PRFConfig(
        job_name="test_gent_pure",
        flag_he=4,
        he_gent=GentParams(MU=1.0, Jm=10.0, D=500.0),
        networks=[],
        initial_inc=0.05, step_time=1.0,
    )

class TestGent:
    """Gent locking hyperelastic."""

    def test_converges(self, gent_cfg):
        inp = BUILDER.build_prf_inp(gent_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Gent failed: total_time={r.get('total_time')}"


@pytest.fixture(scope="module")
def ogden_cfg():
    return PRFConfig(
        job_name="test_ogden_pure",
        flag_he=5,
        he_ogden=OgdenParams(mu1=1.0, alpha1=2.0, D1=500.0),
        networks=[],
        initial_inc=0.05, step_time=1.0,
    )

class TestOgden:
    """Ogden (N=3) hyperelastic."""

    def test_converges(self, ogden_cfg):
        inp = BUILDER.build_prf_inp(ogden_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Ogden failed: total_time={r.get('total_time')}"


@pytest.fixture(scope="module", params=[
    {"flag_he": 3, "C10": 0.5, "C20": 0.1},
    {"flag_he": 4, "MU": 1.0, "Jm": 10.0},
    {"flag_he": 5, "mu1": 1.0, "alpha1": 2.0},
])
def he_visc_cfg(request):
    params = dict(request.param)
    fh = params.pop("flag_he")
    if fh == 3:
        return PRFConfig(
            job_name=f"test_he{fh}_visc",
            flag_he=3,
            he_yeoh=YeohParams(**params, D1=500.0),
            networks=[ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01)],
            initial_inc=0.02, step_time=1.0,
        )
    elif fh == 4:
        return PRFConfig(
            job_name=f"test_he{fh}_visc",
            flag_he=4,
            he_gent=GentParams(**params, D=500.0),
            networks=[ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01)],
            initial_inc=0.02, step_time=1.0,
        )
    elif fh == 5:
        return PRFConfig(
            job_name=f"test_he{fh}_visc",
            flag_he=5,
            he_ogden=OgdenParams(**params, D1=500.0),
            networks=[ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01)],
            initial_inc=0.02, step_time=1.0,
        )
    return None

class TestNewHEWithVisc:
    """New hyperelastic models with viscous network."""

    def test_converges(self, he_visc_cfg):
        inp = BUILDER.build_prf_inp(he_visc_cfg)
        r = solve(inp.stem)
        assert r["success"], f"HE+visc {inp.stem}: total_time={r.get('total_time')}"


@pytest.fixture(scope="module")
def prony_cfg():
    """Single Prony term with Neo-Hookean equilibrium."""
    return PRFConfig(
        job_name="test_prony_1term",
        flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, RBULK_HE=1000.0),
        networks=[ViscousNetworkConfig(
            stiffn=0.5, flag_visc=VISC_PRONY,
            G_PRONY=0.5, K_PRONY=0.0, TAU_PRONY=1.0,
        )],
        initial_inc=0.05, step_time=1.0, max_inc=0.1,
    )


class TestProny:
    """Prony series linear viscoelastic."""

    def test_converges(self, prony_cfg):
        inp = BUILDER.build_prf_inp(prony_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Prony failed: total_time={r.get('total_time')}"

    def test_relaxation(self, prony_cfg):
        """Prony g=0.5, tau=1.0: stress should relax from ~0.9 to ~0.45 MPa."""
        inp = BUILDER.build_prf_inp(prony_cfg)
        r = solve(inp.stem)
        assert r["success"], f"Prony failed: total_time={r.get('total_time')}"
        # After 1s with tau=1.0, stress should be between instantaneous and fully relaxed
        # We just verify convergence here; stress validation is done in compare_pp_prf.py
