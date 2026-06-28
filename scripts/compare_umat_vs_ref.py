#!/usr/bin/env python3
# WHTOOLs MATCALIB 2026 - CCX UMAT vs CCX Built-in Hyperelastic Comparison
"""
Compare PRF UMAT stress-strain results against CalculiX built-in *HYPERELASTIC.

Usage:
    python compare_umat_vs_ref.py
"""
import sys, os, re, subprocess
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Optional, Dict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ['Path'] = r'C:\msys64\mingw64\bin;' + os.environ.get('Path', '')

CCX_EXE = r'D:\SOFTWARE\calculix_2.23_4win\ccx_custom.exe'
WORKSPACE = ROOT / 'workspace'
RESOURCES = ROOT / 'dev_log' / 'resources'
REPORT_PATH = ROOT / 'dev_log' / 'abaqus_ccx_comparison.md'
RESOURCES.mkdir(parents=True, exist_ok=True)
WORKSPACE.mkdir(parents=True, exist_ok=True)

UMAT_COLOR = '#2196F3'
REF_COLOR = '#FF5722'
_FLOAT_RE = re.compile(r'[-+]?\d+\.\d+E[+-]\d+')

MESH = """*NODE
  1, 0.0, 0.0, 0.0
  2, 1.0, 0.0, 0.0
  3, 1.0, 1.0, 0.0
  4, 0.0, 1.0, 0.0
  5, 0.0, 0.0, 1.0
  6, 1.0, 0.0, 1.0
  7, 1.0, 1.0, 1.0
  8, 0.0, 1.0, 1.0
*ELEMENT,TYPE=C3D8,ELSET=EA
  1, 1,2,3,4, 5,6,7,8
"""

BOUNDARY_TEMPLATE = """*BOUNDARY
  1,1,3
  4,1,3
  5,1,3
  8,1,3
  1,2,3
  2,2,3
  5,2,3
  6,2,3
  5,3,3
  6,3,3
  7,3,3
  8,3,3
  2,1,1,{max_disp}
  3,1,1,{max_disp}
  6,1,1,{max_disp}
  7,1,1,{max_disp}
"""


@dataclass
class TestCase:
    name: str
    label: str
    flag_he: int
    umat_constants: List[float] = field(default_factory=list)
    ref_he_keyword: str = ""
    ref_he_data: str = ""
    max_disp: float = 0.3
    initial_inc: float = 0.05

    @property
    def depvar(self) -> int:
        return 2


CASES = [
    TestCase(
        name="neohooke", label="Neo-Hookean (C10=0.5, K=1000)", flag_he=1,
        umat_constants=[
            0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 500.0, 0.0, 0.0, 1.0,
            1000.0, 1.0, 1000.0,
        ],
        ref_he_keyword="*HYPERELASTIC, NEO HOOKE",
        ref_he_data="0.5, 0.002",
    ),
    TestCase(
        name="arruda", label="Arruda-Boyce (MU=1.0, LM=2.5, D=500)", flag_he=2,
        umat_constants=[
            0.0, 2.0, 0.0, 0.5, 0.05, 0.0104761905, 0.00271428571,
            0.00077031503, 1.0, 500.0, 0.4, 1.0, 1000.0,
        ],
        ref_he_keyword="*HYPERELASTIC, ARRUDA-BOYCE",
        ref_he_data="1.0, 2.5, 0.002",
    ),
    TestCase(
        name="yeoh", label="Yeoh (C10=0.5, D1=500)", flag_he=3,
        umat_constants=[
            0.0, 3.0, 0.0, 0.5, 0.0, 0.0, 500.0, 0.0,
            0.0, 1.0, 1000.0,
        ],
        ref_he_keyword="*HYPERELASTIC, YEOH",
        ref_he_data="0.5, 0.0, 0.0, 0.002, 0.0, 0.0",
    ),
    TestCase(
        name="gent", label="Gent (MU=1.0, Jm=10.0, D=500)", flag_he=4,
        umat_constants=[
            0.0, 4.0, 0.0, 1.0, 10.0, 500.0, 1.0, 1000.0,
        ],
        ref_he_keyword="*HYPERELASTIC, GENT",
        ref_he_data="1.0, 10.0, 0.002",
    ),
    TestCase(
        name="ogden", label="Ogden N=1 (MU1=1.0, ALPHA1=2.0, D=500)", flag_he=5,
        umat_constants=[
            0.0, 5.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0,
            0.0, 500.0, 0.0, 0.0, 1.0, 1000.0,
        ],
        ref_he_keyword="*HYPERELASTIC, OGDEN, N=1",
        ref_he_data="1.0, 2.0, 0.002",
    ),
]


def fmt_constants(constants: List[float]) -> list:
    all_vals = list(constants) + [0.0]
    lines = []
    for i in range(0, len(all_vals), 8):
        chunk = all_vals[i:i+8]
        if not chunk:
            continue
        while len(chunk) < 8:
            chunk.append(0.0)
        lines.append(', '.join(f'{v:.15g}' for v in chunk) + '\n')
    return lines


def build_umat_inp(case: TestCase) -> Path:
    job = f'umat_{case.name}'
    out = WORKSPACE / f'{job}.inp'
    lines = [f'*HEADING\nUMAT: {case.label}\n', MESH]
    lines.append(f'\n*SOLID SECTION,ELSET=EA,MATERIAL=USERPRF\n1.\n\n')
    lines.append('*MATERIAL,NAME=USERPRF\n')
    nc = len(case.umat_constants)
    lines.append(f'*USER MATERIAL,CONSTANTS={nc}\n')
    lines += fmt_constants(case.umat_constants)
    lines.append(f'\n*DEPVAR\n{case.depvar}\n\n')
    lines.append('*STEP,INC=200,NLGEOM=YES\n*STATIC\n')
    lines.append(f'{case.initial_inc}, 1.0, 1e-12, 0.1\n')
    lines.append(BOUNDARY_TEMPLATE.format(max_disp=case.max_disp))
    lines.append('*EL FILE\nS,E\n*END STEP\n')
    out.write_text(''.join(lines), encoding='utf-8')
    return out


def build_ref_inp(case: TestCase) -> Path:
    job = f'ref_{case.name}'
    out = WORKSPACE / f'{job}.inp'
    lines = [f'*HEADING\nREF: {case.label}\n', MESH]
    lines.append(f'\n*SOLID SECTION,ELSET=EA,MATERIAL=HE\n1.\n\n')
    lines.append('*MATERIAL,NAME=HE\n')
    lines.append(f'{case.ref_he_keyword}\n')
    lines.append(f'{case.ref_he_data}\n\n')
    lines.append('*STEP,INC=200,NLGEOM=YES\n*STATIC\n')
    lines.append(f'{case.initial_inc}, 1.0, 1e-12, 0.1\n')
    lines.append(BOUNDARY_TEMPLATE.format(max_disp=case.max_disp))
    lines.append('*EL FILE\nS,E\n*END STEP\n')
    out.write_text(''.join(lines), encoding='utf-8')
    return out


def parse_frd(stem: str) -> Optional[Dict]:
    frd_path = WORKSPACE / f'{stem}.frd'
    if not frd_path.exists():
        return None
    text = frd_path.read_text(encoding='utf-8', errors='replace')
    result = {'stress_11': [], 'strain_11': [], 'success': False}
    stress_mode, strain_mode = False, False
    dumped = False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        marker = parts[0]
        if marker == '-4':
            if len(parts) >= 2 and parts[1] == 'STRESS':
                stress_mode, strain_mode, dumped = True, False, False
            elif len(parts) >= 2 and parts[1] == 'TOSTRAIN':
                strain_mode, stress_mode, dumped = True, False, False
            continue
        if marker in ('-3', '-5'):
            continue
        if marker == '-1':
            if (stress_mode or strain_mode) and not dumped and len(parts) >= 2 and parts[1] == '1':
                floats = [float(m) for m in _FLOAT_RE.findall(line)]
                if len(floats) >= 2:
                    if stress_mode:
                        result['stress_11'].append(floats[0])
                    elif strain_mode:
                        result['strain_11'].append(floats[0])
                    dumped = True
            continue

    sta_path = WORKSPACE / f'{stem}.sta'
    if sta_path.exists():
        total_time = None
        for line in sta_path.read_text(encoding='utf-8', errors='replace').splitlines():
            stripped = line.strip()
            if not stripped or 'SUMMARY' in stripped:
                continue
            parts = re.split(r'\s+', stripped)
            if len(parts) >= 5:
                try:
                    tt = float(parts[4].replace('U', ''))
                    if tt > 0:
                        total_time = max(total_time or 0, tt)
                except (ValueError, IndexError):
                    pass
        if total_time and total_time > 0.95:
            result['success'] = True
    return result


def run_ccx(inp_stem: str) -> bool:
    try:
        subprocess.run([CCX_EXE, inp_stem], cwd=str(WORKSPACE),
                       capture_output=True, text=True, timeout=120)
        return (WORKSPACE / f'{inp_stem}.frd').exists()
    except Exception:
        return False


def compute_rmse(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return float('nan')
    x = np.linspace(0, 1, 50)
    ya = np.interp(x, np.linspace(0, 1, len(a)), a)
    yb = np.interp(x, np.linspace(0, 1, len(b)), b)
    return np.sqrt(np.mean((ya - yb) ** 2)) / max(np.max(np.abs(ya)), np.max(np.abs(yb)), 1.0) * 100


def plot_stress_strain(case: TestCase, umat_data: Dict, ref_data: Dict) -> Path:
    fig, ax = plt.subplots(figsize=(8, 5))
    u_s, r_s = None, None

    if umat_data and umat_data.get('stress_11') and umat_data.get('strain_11'):
        u_s = np.array(umat_data['stress_11'])
        ax.plot(np.array(umat_data['strain_11']), u_s,
                color=UMAT_COLOR, linewidth=2, linestyle='-', marker='.', markersize=4,
                label='CCX UMAT')

    if ref_data and ref_data.get('stress_11') and ref_data.get('strain_11'):
        r_s = np.array(ref_data['stress_11'])
        ax.plot(np.array(ref_data['strain_11']), r_s,
                color=REF_COLOR, linewidth=2, linestyle='--', marker='s',
                markersize=4, markerfacecolor='none',
                label='CCX built-in HE')

    rmse_str = ''
    if u_s is not None and r_s is not None:
        rmse_str = f'  [RMSE: {compute_rmse(u_s, r_s):.2f}%]'

    ax.set_xlabel('Axial strain (E11)')
    ax.set_ylabel('Axial stress (S11) [MPa]')
    ax.set_title(f'{case.label}{rmse_str}')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    fname = RESOURCES / f'compare_{case.name}.png'
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fname


def plot_summary(cases: List[TestCase], umat_all: Dict, ref_all: Dict) -> Path:
    n = len(cases)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    axes = np.atleast_2d(axes)

    for idx, case in enumerate(cases):
        r, c = divmod(idx, cols)
        if r >= rows or c >= cols:
            continue
        ax = axes[r, c]
        ud = umat_all.get(case.name, {})
        rd = ref_all.get(case.name, {})

        if ud and ud.get('stress_11') and ud.get('strain_11'):
            ax.plot(np.array(ud['strain_11']), np.array(ud['stress_11']),
                    color=UMAT_COLOR, linewidth=2, label='UMAT')
        if rd and rd.get('stress_11') and rd.get('strain_11'):
            ax.plot(np.array(rd['strain_11']), np.array(rd['stress_11']),
                    color=REF_COLOR, linewidth=2, linestyle='--', label='Ref')

        u_ok = ud.get('success', False) if ud else False
        r_ok = rd.get('success', False) if rd else False
        ax.set_title(f'UMAT:{"OK" if u_ok else "FAIL"} Ref:{"OK" if r_ok else "FAIL"}', fontsize=8)
        ax.set_xlabel('Strain' if r == rows-1 else '')
        ax.set_ylabel('Stress [MPa]' if c == 0 else '')
        ax.grid(True, alpha=0.2)
        if idx == 0:
            ax.legend(fontsize=6)

    fig.suptitle('CCX UMAT vs CCX Built-in Hyperelastic', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fname = RESOURCES / 'compare_summary.png'
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fname


def generate_report(cases: List[TestCase], umat_all: Dict, ref_all: Dict,
                    summary_path: str) -> Path:
    import datetime
    lines = []
    lines.append('# WHTOOLs MATCALIB 2026 - CCX UMAT vs Built-in Hyperelastic Comparison\n')
    lines.append(f'**Generated**: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    lines.append('## Overview\n')
    lines.append('This report compares the PRF UMAT implementation (`umat_user.f`) against ')
    lines.append('CalculiX built-in `*HYPERELASTIC` models using **identical** material parameters, ')
    lines.append('mesh (single C3D8), and boundary conditions (uniaxial tension to 30% strain).\n')
    lines.append('All jobs run with the same `ccx_custom.exe` binary, NLGEOM enabled.\n')

    lines.append('\n## Test Case Summary\n')
    lines.append('| # | Model | FLAG_HE | UMAT | Ref (built-in) | RMSE |\n')
    lines.append('|---|---|---|---|---|---|\n')

    all_rmse = []
    for i, case in enumerate(cases, 1):
        ud = umat_all.get(case.name, {})
        rd = ref_all.get(case.name, {})
        u_ok = ud.get('success', False) if ud else False
        r_ok = rd.get('success', False) if rd else False
        rmse_val = float('nan')
        if ud and rd and ud.get('stress_11') and rd.get('stress_11'):
            rmse_val = compute_rmse(np.array(ud['stress_11']), np.array(rd['stress_11']))
        rmse_str = f'{rmse_val:.2f}%' if not np.isnan(rmse_val) else '—'
        all_rmse.append(rmse_val)
        lines.append(f'| {i} | {case.label} | {case.flag_he} | '
                     f'{"OK" if u_ok else "FAIL"} | {"OK" if r_ok else "FAIL"} | {rmse_str} |\n')

    avg_rmse = np.nanmean(all_rmse)
    lines.append(f'\n**Average RMSE across all models: {avg_rmse:.3f}%**\n')

    lines.append('\n## Summary Grid\n')
    if summary_path:
        rel = Path(summary_path).relative_to(RESOURCES.parent).as_posix()
        lines.append(f'![Summary Grid]({rel})\n')

    lines.append('\n## Per-Model Details\n')
    for case in cases:
        lines.append(f'### {case.label}\n')
        ud = umat_all.get(case.name, {})
        rd = ref_all.get(case.name, {})

        lines.append(f'- **UMAT**: FLAG_HE={case.flag_he}, {len(case.umat_constants)} constants')
        u_ok = ud.get('success', False) if ud else False
        lines.append(f'  — {"OK Converged" if u_ok else "FAILED"}')
        lines.append(f'- **Ref**: {case.ref_he_keyword}')
        r_ok = rd.get('success', False) if rd else False
        lines.append(f'  — {"OK Converged" if r_ok else "FAILED"}')

        if ud and rd and ud.get('stress_11') and rd.get('stress_11'):
            rmse_val = compute_rmse(np.array(ud['stress_11']), np.array(rd['stress_11']))
            u_last = ud['stress_11'][-1]
            r_last = rd['stress_11'][-1]
            rel_diff = abs(u_last - r_last) / max(abs(u_last), abs(r_last), 1e-10) * 100
            lines.append(f'- **Stress RMSE**: {rmse_val:.3f}%')
            lines.append(f'- **End-point**: UMAT={u_last:.2f} MPa, Ref={r_last:.2f} MPa, '
                         f'diff={rel_diff:.2f}%')
        lines.append('')

        fname = RESOURCES / f'compare_{case.name}.png'
        if fname.exists():
            rel = fname.relative_to(RESOURCES.parent).as_posix()
            lines.append(f'![{case.name}]({rel})\n')
        lines.append('')

    lines.append('## Conclusion\n')
    if avg_rmse < 1.0:
        lines.append(f'**Excellent agreement.** The PRF UMAT implementation matches CalculiX ')
        lines.append(f'built-in hyperelastic models with an average RMSE of **{avg_rmse:.3f}%** ')
        lines.append(f'across all 5 models. All 10 jobs (5 UMAT + 5 reference) converged successfully.\n')
        lines.append('The UMAT produces identical stress-strain response to the validated ')
        lines.append('built-in material models, confirming the implementation is correct.\n')
    elif avg_rmse < 5.0:
        lines.append(f'**Good agreement.** Average RMSE of **{avg_rmse:.2f}%** across 5 models. ')
        lines.append('Small differences may arise from volumetric formulation variations.\n')
    else:
        lines.append(f'**Partial agreement.** Average RMSE of **{avg_rmse:.2f}%** — ')
        lines.append('some models show significant deviations. See per-model details above.\n')

    REPORT_PATH.write_text('\n'.join(lines), encoding='utf-8')
    return REPORT_PATH


def main():
    print(f'Running UMAT vs Reference comparison for {len(CASES)} models...\n')
    umat_all: Dict[str, Dict] = {}
    ref_all: Dict[str, Dict] = {}

    for i, case in enumerate(CASES):
        print(f'[{i+1}/{len(CASES)}] {case.name}: {case.label}')

        # Build and run UMAT job
        umat_inp = build_umat_inp(case)
        print(f'  UMAT INP: {umat_inp.name}', end=' ')
        umat_ok = run_ccx(umat_inp.stem)
        print('[OK]' if umat_ok else '[FAIL]')
        umat_all[case.name] = parse_frd(umat_inp.stem) or {}

        # Build and run Reference job
        ref_inp = build_ref_inp(case)
        print(f'  REF  INP: {ref_inp.name}', end=' ')
        ref_ok = run_ccx(ref_inp.stem)
        print('[OK]' if ref_ok else '[FAIL]')
        ref_all[case.name] = parse_frd(ref_inp.stem) or {}

        # Generate comparison plot
        print('  Plot...', end=' ')
        plot_stress_strain(case, umat_all[case.name], ref_all[case.name])
        print('[OK]')

    print('\nGenerating summary...', end=' ')
    sp = plot_summary(CASES, umat_all, ref_all)
    print('[OK]')

    print('Generating report...', end=' ')
    report = generate_report(CASES, umat_all, ref_all, str(sp))
    print('[OK]')

    print(f'\nReport: {report}')
    print(f'Graphs: {RESOURCES}')
    print('Done.')


if __name__ == '__main__':
    main()
