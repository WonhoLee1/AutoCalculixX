#!/usr/bin/env python3
# WHTOOLs MATCALIB 2026 — PRF UMAT stress-strain curve plotter
"""Plot PRF UMAT stress-strain curves from CalculiX .frd output."""
import sys, subprocess, os, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.core.config import (
    PRFConfig, ViscousNetworkConfig, PolynomialHyperelasticParams,
    ArrudaBoyceParams,
)
from src.core.model_builder import CalculixModelBuilder

os.environ['Path'] = r'C:\msys64\mingw64\bin;' + os.environ.get('Path', '')
CCX = r'D:\SOFTWARE\calculix_2.23_4win\ccx_custom.exe'
WORKSPACE = Path(__file__).resolve().parent.parent / 'workspace'


def build_single_step_inp(cfg: PRFConfig) -> Path:
    builder = CalculixModelBuilder(WORKSPACE)
    coords = [
        "0.0, 0.0, 0.0", "1.0, 0.0, 0.0", "1.0, 1.0, 0.0", "0.0, 1.0, 0.0",
        "0.0, 0.0, 1.0", "1.0, 0.0, 1.0", "1.0, 1.0, 1.0", "0.0, 1.0, 1.0",
    ]
    inp = WORKSPACE / f"{cfg.job_name}.inp"
    max_disp = cfg.max_disp
    lines = [f"*HEADING\nPRF: {cfg.job_name}\n", "*NODE\n"]
    for i, xyz in enumerate(coords, 1):
        lines.append(f"{i}, {xyz}\n")
    nids = list(range(1, 9))
    lines.append(f"*ELEMENT,TYPE={cfg.element_type},ELSET=EL1\n")
    lines.append(f"1, {', '.join(map(str, nids))}\n\n")
    lines.append(f"*SOLID SECTION,ELSET=EL1,MATERIAL={cfg.material_name}\n\n")
    lines.append(f"*MATERIAL,NAME={cfg.material_name}\n")
    nc = cfg.nconstants
    lines.append(f"*USER MATERIAL,CONSTANTS={nc}\n")
    cons = list(cfg.constants)
    cons.append(0.0)
    for i in range(0, nc + 1, 8):
        chunk = cons[i:i+8]
        if not chunk:
            continue
        while len(chunk) < 8:
            chunk.append(0.0)
        lines.append(", ".join(f"{v:.15g}" for v in chunk) + "\n")
    lines.append("\n")
    lines.append(f"*DEPVAR\n{cfg.depvar}\n\n")
    lines.append(f"*STEP,INC=500,NLGEOM=YES\n")
    lines.append("*STATIC\n")
    lines.append(f"{cfg.initial_inc}, 1.0, {cfg.min_inc}, {cfg.max_inc}\n")
    lines.append("*BOUNDARY\n")
    lines.append("1,1,3\n2,2,3\n4,2,3\n5,1,3\n6,2,3\n8,2,3\n")
    lines.append(f"1,1,1\n2,1,1,{max_disp}\n3,1,1,{max_disp}\n4,1,1\n")
    lines.append(f"5,1,1\n6,1,1,{max_disp}\n7,1,1,{max_disp}\n8,1,1\n")
    lines.append("*EL FILE\nS,E\n")
    lines.append("*END STEP\n")
    with open(inp, 'w', encoding='utf-8') as f:
        f.writelines(lines)
    return inp


# regex for Fortran-format floats: optional sign, digits, decimal, E+/-exp
_FLOAT_RE = re.compile(r'[-+]?\d+\.\d+E[+-]\d+')

def parse_frd(frd_path: Path):
    if not frd_path.exists():
        return [], []
    text = frd_path.read_text(encoding='utf-8', errors='replace')
    strain, stress = [], []
    lines = text.splitlines()
    in_stress = False
    in_strain = False
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        marker = parts[0]
        if marker == '-4':
            if parts[1] == 'STRESS':
                in_stress = True
                in_strain = False
            elif parts[1] == 'TOSTRAIN':
                in_strain = True
                in_stress = False
            continue
        if marker == '-3':
            in_stress = False
            in_strain = False
            continue
        if marker == '-5':
            continue
        if marker == '-1':
            # Extract all float values using regex (handles fused Fortran output)
            floats = [float(m) for m in _FLOAT_RE.findall(line)]
            if len(floats) < 6:
                continue
            if in_stress and parts[1] == '1':
                stress.append(floats[0])
            elif in_strain and parts[1] == '1':
                strain.append(floats[0])
    return stress, strain


configs = [
    ("N0 (Pure HE) C10=0.5", PRFConfig(
        job_name="plot_n0", flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[], initial_inc=0.05, step_time=1.0, max_inc=0.1, max_disp=0.3)),
    ("1 BB A1=0.1 EXPC=0", PRFConfig(
        job_name="plot_1bb", flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[ViscousNetworkConfig(stiffn=0.5, A1=0.1, EXPC=0.0, EXPM=2.0, KSI=0.01)],
        initial_inc=0.02, step_time=1.0, max_inc=0.05, max_disp=0.15)),
    ("1 BB A1=0.01 EXPC=-0.5", PRFConfig(
        job_name="plot_1bb_slow", flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=-0.5, EXPM=3.0, KSI=0.001)],
        initial_inc=0.02, step_time=1.0, max_inc=0.1, max_disp=0.2)),
    ("1 BB A1=0.5 EXPC=1", PRFConfig(
        job_name="plot_1bb_fast", flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[ViscousNetworkConfig(stiffn=0.8, A1=0.5, EXPC=1.0, EXPM=1.5, KSI=0.01)],
        initial_inc=0.02, step_time=1.0, max_inc=0.1, max_disp=0.2)),
    ("2 BB A1=0.01+0.005", PRFConfig(
        job_name="plot_2bb", flag_he=1,
        he_poly=PolynomialHyperelasticParams(C10=0.5, C01=0.0, RBULK_HE=1000.0),
        networks=[
            ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01),
            ViscousNetworkConfig(stiffn=0.2, A1=0.005, EXPC=0.0, EXPM=2.0, KSI=0.01),
        ],
        initial_inc=0.01, step_time=1.0, max_inc=0.1, max_disp=0.2)),
    ("Arruda-Boyce + BB", PRFConfig(
        job_name="plot_ab_bb", flag_he=2,
        he_ab=ArrudaBoyceParams(C_R=0.5, K0=500.0, N=2.5),
        networks=[ViscousNetworkConfig(stiffn=0.3, A1=0.01, EXPC=0.0, EXPM=2.0, KSI=0.01)],
        initial_inc=0.01, step_time=1.0, max_inc=0.1, max_disp=0.2)),
]

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
print(f"{'Config':30s} {'stress':>6s} {'strain':>6s} {'S11_last':>10s} {'E11_last':>10s}")
print("-"*70)

for idx, (label, cfg) in enumerate(configs):
    ax = axes[idx]
    inp = build_single_step_inp(cfg)
    for f in [inp.with_suffix('.sta'), inp.with_suffix('.frd')]:
        if f.exists(): f.unlink()

    result = subprocess.run([CCX, inp.stem], cwd=str(WORKSPACE),
                            capture_output=True, text=True, timeout=300)

    stress, strain = parse_frd(inp.with_suffix('.frd'))
    n_pts = len(stress)

    # only plot if we have 3+ points
    if n_pts >= 3:
        ax.plot(strain, stress, 'o-', color=colors[idx % len(colors)],
                linewidth=2, markersize=4, label=label)
        ax.set_xlabel('E11 (Log Strain)', fontsize=12)
        ax.set_ylabel('S11 (Cauchy Stress)', fontsize=12)
        ax.set_title(f"{label}\n({n_pts} inc)", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
    else:
        ax.text(0.5, 0.5, f"NO DATA\n({n_pts} pts)", ha='center', va='center',
                transform=ax.transAxes, fontsize=12)
        ax.set_title(label, fontsize=10)

    last_s = stress[-1] if stress else -1
    last_e = strain[-1] if strain else -1
    print(f"{label:30s} {n_pts:6d} {len(stress):6d} {last_s:10.4f} {last_e:10.6f}")

plt.tight_layout()
out = WORKSPACE / 'prf_stress_strain_curves.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out}")
plt.close()
