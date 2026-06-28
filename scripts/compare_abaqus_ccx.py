#!/usr/bin/env python3
# WHTOOLs MATCALIB 2026 - Abaqus vs CalculiX PRF UMAT Comparison
"""
Compare CalculiX UMAT stress-strain against Abaqus built-in *HYPERELASTIC.

D convention: CCX UMAT uses dU/dJ = 2*D*(J-1) (K₀ = 2D, D=500 -> K=1000)
              Abaqus uses  dU/dJ = (2/D)*(J-1) (K₀ = 2/D, D=0.002 -> K=1000)
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

CCX_COLOR = '#2196F3'
ABQ_COLOR = '#FF5722'
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

BOUNDARY = """*BOUNDARY
1,1,3
4,1,1
4,3,3
8,1,1
5,1,2
2,2,3
6,2,2
2,3,3
3,3,3
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
    ccx_umat_constants: List[float]
    abq_he_keyword: str
    abq_he_data: str
    max_disp: float = 0.3
    initial_inc: float = 0.05


CASES = [
    TestCase(
        name="neohooke", label="Neo-Hookean", flag_he=1,
        ccx_umat_constants=[
            0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 500.0, 0.0, 0.0, 1.0,
            1000.0, 1.0, 1000.0,
        ],
        abq_he_keyword="*HYPERELASTIC,NEO HOOKE",
        abq_he_data="0.5, 0.002",
    ),
    TestCase(
        name="arruda", label="Arruda-Boyce", flag_he=2,
        ccx_umat_constants=[
            0.0, 2.0, 0.0, 0.5, 0.05, 0.0104761905, 0.00271428571,
            0.00077031503, 1.0, 500.0, 0.16, 1.0, 1000.0,
        ],
        abq_he_keyword="*HYPERELASTIC,ARRUDA-BOYCE",
        abq_he_data="1.0, 2.5, 0.002",
    ),
    TestCase(
        name="yeoh", label="Yeoh", flag_he=3,
        ccx_umat_constants=[
            0.0, 3.0, 0.0, 0.5, 0.0, 0.0, 500.0, 0.0,
            0.0, 1.0, 1000.0,
        ],
        abq_he_keyword="*HYPERELASTIC,YEOH",
        abq_he_data="0.5, 0.0, 0.0, 0.002, 0.0, 0.0",
    ),
    TestCase(
        name="ogden", label="Ogden N=1", flag_he=5,
        ccx_umat_constants=[
            0.0, 5.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0,
            0.0, 500.0, 0.0, 0.0, 1.0, 1000.0,
        ],
        abq_he_keyword="*HYPERELASTIC,OGDEN,N=1",
        abq_he_data="1.0, 2.0, 0.002",
    ),
]


def fmt_8(constants: List[float]) -> list:
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


# ─── CCX UMAT INP builder ──────────────────────────────────────────────────

def build_ccx_inp(case: TestCase) -> Path:
    job = f'ccx_{case.name}'
    out = WORKSPACE / f'{job}.inp'
    lines = [f'*HEADING\nCCX UMAT: {case.label}\n', MESH,
             f'*SOLID SECTION,ELSET=EA,MATERIAL=USERPRF\n1.\n\n',
             '*MATERIAL,NAME=USERPRF\n']
    nc = len(case.ccx_umat_constants)
    lines.append(f'*USER MATERIAL,CONSTANTS={nc}\n')
    lines += fmt_8(case.ccx_umat_constants)
    lines.append('\n*DEPVAR\n2\n\n')
    lines.append('*STEP,INC=200,NLGEOM=YES\n*STATIC\n')
    lines.append(f'{case.initial_inc}, 1.0, 1e-12, 0.1\n')
    lines.append(BOUNDARY.format(max_disp=case.max_disp))
    lines.append('*EL FILE\nS,E\n*END STEP\n')
    out.write_text(''.join(lines), encoding='utf-8')
    return out


# ─── Abaqus INP builder ────────────────────────────────────────────────────

def build_abq_inp(case: TestCase) -> Path:
    job = f'abq_{case.name}'
    out = WORKSPACE / f'{job}.inp'
    lines = [f'*HEADING\nAbaqus: {case.label}\n', MESH,
             f'*SOLID SECTION,ELSET=EA,MATERIAL=HE\n1.\n\n',
             '*MATERIAL,NAME=HE\n',
             f'{case.abq_he_keyword}\n',
             f'{case.abq_he_data}\n',
             '*DIAGNOSTICS,NONHYBRID=WARNING\n\n',
             '*STEP,INC=200,NLGEOM=YES\n*STATIC\n',
             f'{case.initial_inc}, 1.0, 1e-12, 0.1\n',
             BOUNDARY.format(max_disp=case.max_disp),
             '*EL PRINT,ELSET=EA,FREQ=1,POSITION=CENTROIDAL\nS,\n',
             '*END STEP\n']
    out.write_text(''.join(lines), encoding='utf-8')
    return out


# ─── Parser: CCX .frd ─────────────────────────────────────────────────────

def parse_ccx(stem: str) -> Optional[Dict]:
    frd = WORKSPACE / f'{stem}.frd'
    if not frd.exists():
        return None
    text = frd.read_text(encoding='utf-8', errors='replace')
    result = {'stress_11': [], 'strain_11': [], 'success': False}
    smode, emode, dumped = False, False, False
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        p = s.split()
        if p[0] == '-4':
            if len(p) >= 2 and p[1] == 'STRESS':
                smode, emode, dumped = True, False, False
            elif len(p) >= 2 and p[1] == 'TOSTRAIN':
                emode, smode, dumped = True, False, False
            continue
        if p[0] in ('-3', '-5'):
            continue
        if p[0] == '-1' and (smode or emode) and not dumped and len(p) >= 2 and p[1] == '1':
            floats = [float(m) for m in _FLOAT_RE.findall(line)]
            if len(floats) >= 2:
                (result['stress_11'] if smode else result['strain_11']).append(floats[0])
                dumped = True
        continue
    sta = WORKSPACE / f'{stem}.sta'
    if sta.exists():
        tt = None
        for line in sta.read_text(encoding='utf-8', errors='replace').splitlines():
            s = line.strip()
            if not s or 'SUMMARY' in s:
                continue
            parts = re.split(r'\s+', s)
            if len(parts) >= 5:
                try:
                    t = float(parts[4].replace('U', ''))
                    if t > 0:
                        tt = max(tt or 0, t)
                except (ValueError, IndexError):
                    pass
        if tt and tt > 0.95:
            result['success'] = True
    return result


# ─── Parser: Abaqus .dat ──────────────────────────────────────────────────

def parse_abq(stem: str) -> Optional[Dict]:
    dat = WORKSPACE / f'{stem}.dat'
    if not dat.exists():
        return None
    text = dat.read_text(encoding='utf-8', errors='replace')
    result = {'stress_11': [], 'strain_11': [], 'success': False}
    stress_mode = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith('ELEMENT') and 'S11' in s:
            stress_mode = True
            continue
        if stress_mode and (s.startswith('MAXIMUM') or s.startswith('MINIMUM')
                           or s.startswith('THE TABLE')):
            stress_mode = False
            continue
        if stress_mode:
            parts = s.split()
            if len(parts) >= 3 and parts[0] == '1':
                try:
                    result['stress_11'].append(float(parts[1].replace('D', 'E')))
                except (ValueError, IndexError):
                    pass

    sta = WORKSPACE / f'{stem}.sta'
    if sta.exists():
        for line in sta.read_text(encoding='utf-8', errors='replace').splitlines():
            if 'COMPLETED SUCCESSFULLY' in line:
                result['success'] = True
                break
    return result


# ─── Runners ──────────────────────────────────────────────────────────────

def run_ccx(stem: str) -> bool:
    try:
        subprocess.run([CCX_EXE, stem], cwd=str(WORKSPACE),
                       capture_output=True, text=True, timeout=120)
        return (WORKSPACE / f'{stem}.frd').exists()
    except Exception:
        return False


def run_abq(stem: str) -> bool:
    try:
        subprocess.run(['abaqus', f'job={stem}', f'input={stem}.inp', 'interactive'],
                       cwd=str(WORKSPACE), timeout=600, shell=True)
        return (WORKSPACE / f'{stem}.odb').exists()
    except Exception:
        return False


# ─── Plotting ─────────────────────────────────────────────────────────────

def compute_rmse(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return float('nan')
    x = np.linspace(0, 1, 50)
    ya = np.interp(x, np.linspace(0, 1, len(a)), a)
    yb = np.interp(x, np.linspace(0, 1, len(b)), b)
    return np.sqrt(np.mean((ya - yb)**2)) / max(np.max(np.abs(ya)), np.max(np.abs(yb)), 1.0) * 100


def plot(case: TestCase, ccx_d: Dict, abq_d: Dict) -> Path:
    fig, ax = plt.subplots(figsize=(9, 5))
    cs, a_s = None, None

    if ccx_d and ccx_d.get('stress_11') and ccx_d.get('strain_11'):
        cs = np.array(ccx_d['stress_11'])
        ax.plot(np.array(ccx_d['strain_11']), cs, color=CCX_COLOR, linewidth=2,
                linestyle='-', marker='.', markersize=4, label='CalculiX UMAT')

    if abq_d and abq_d.get('stress_11'):
        a_s = np.array(abq_d['stress_11'])
        n = len(a_s)
        abq_strain = np.linspace(0, case.max_disp, n)
        ax.plot(abq_strain, a_s, color=ABQ_COLOR, linewidth=2, linestyle='--',
                marker='s', markersize=4, markerfacecolor='none',
                label='Abaqus built-in HE')

    rmse_str = ''
    if cs is not None and a_s is not None:
        rmse_str = f'  [RMSE: {compute_rmse(cs, a_s):.2f}%]'

    ax.set_xlabel('Axial strain (E11)')
    ax.set_ylabel('Axial stress (S11) [MPa]')
    ax.set_title(f'{case.label} — CCX UMAT vs Abaqus{rmse_str}')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    fname = RESOURCES / f'abq_ccx_{case.name}.png'
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fname


def plot_summary(cases: List[TestCase], ccx_all: Dict, abq_all: Dict) -> Path:
    n = len(cases)
    cols = min(2, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4.5*rows))
    axes = np.atleast_2d(axes)

    for idx, case in enumerate(cases):
        r, c = divmod(idx, cols)
        if r >= rows or c >= cols:
            continue
        ax = axes[r, c]
        cd = ccx_all.get(case.name, {})
        ad = abq_all.get(case.name, {})

        if cd and cd.get('stress_11') and cd.get('strain_11'):
            ax.plot(np.array(cd['strain_11']), np.array(cd['stress_11']),
                    color=CCX_COLOR, linewidth=2, label='CCX')
        if ad and ad.get('stress_11'):
            a_s = np.array(ad['stress_11'])
            ax.plot(np.linspace(0, case.max_disp, len(a_s)), a_s,
                    color=ABQ_COLOR, linewidth=2, linestyle='--', label='ABQ')

        c_ok = cd.get('success', False) if cd else False
        a_ok = ad.get('success', False) if ad else False
        rmse_val = float('nan')
        if cd and ad and cd.get('stress_11') and ad.get('stress_11'):
            rmse_val = compute_rmse(np.array(cd['stress_11']), np.array(ad['stress_11']))
        rmse_s = f'RMSE:{rmse_val:.2f}%' if not np.isnan(rmse_val) else ''
        ax.set_title(f'CCX:{"OK" if c_ok else "FAIL"} ABQ:{"OK" if a_ok else "FAIL"}  {rmse_s}', fontsize=8)
        ax.set_xlabel('Strain' if r == rows-1 else '')
        ax.set_ylabel('Stress [MPa]' if c == 0 else '')
        ax.grid(True, alpha=0.2)
        if idx == 0:
            ax.legend(fontsize=6)

    fig.suptitle('CalculiX UMAT vs Abaqus Built-in Hyperelastic\nSingle C3D8, Uniaxial Tension 30%', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fname = RESOURCES / 'abq_ccx_summary.png'
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fname


# ─── Report ───────────────────────────────────────────────────────────────

def generate_report(cases: List[TestCase], ccx_all: Dict, abq_all: Dict, summary: str) -> Path:
    import datetime
    lines = [
        '# WHTOOLs MATCALIB 2026 — Abaqus vs CalculiX UMAT Comparison\n',
        f'**Generated**: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n',
        '## Overview\n',
        'This report compares the PRF UMAT (`umat_user.f`) results from CalculiX against Abaqus built-in `*HYPERELASTIC`.\n',
        '**Identical**: single C3D8 element, uniaxial tension to 30% strain, NLGEOM.\n',
        '**Material**: same C10=0.5, K=1000 (CCX D=500 → ABQ D=0.002).\n',
        '> **Note**: Abaqus uses `*DIAGNOSTICS,NONHYBRID=WARNING` for nearly-incompressible C3D8 (ν≈0.4995).\n',
        '\n## Results Summary\n',
        '| # | Model | FLAG_HE | CalculiX UMAT | Abaqus | RMSE | End-point diff |\n',
        '|---|---|---|---|---|---|---|\n',
    ]

    all_rmse = []
    for i, case in enumerate(cases, 1):
        cd = ccx_all.get(case.name, {})
        ad = abq_all.get(case.name, {})
        c_ok = cd.get('success', False) if cd else False
        a_ok = ad.get('success', False) if ad else False
        rmse_val = float('nan')
        ep_diff = float('nan')
        if cd and ad and cd.get('stress_11') and ad.get('stress_11'):
            cs = np.array(cd['stress_11'])
            a_s = np.array(ad['stress_11'])
            rmse_val = compute_rmse(cs, a_s)
            ep_diff = abs(cs[-1] - a_s[-1]) / max(abs(cs[-1]), abs(a_s[-1]), 1e-10) * 100
        all_rmse.append(rmse_val)
        lines.append(f'| {i} | {case.label} | {case.flag_he} | {"OK" if c_ok else "FAIL"} | '
                     f'{"OK" if a_ok else "FAIL"} | '
                     f'{rmse_val:.2f}%' if not np.isnan(rmse_val) else '—')
        if not np.isnan(ep_diff):
            lines[-1] += f' | {ep_diff:.2f}% |\n'
        else:
            lines[-1] += ' | — |\n'

    avg = np.nanmean(all_rmse)
    lines.append(f'\n**Average RMSE: {avg:.3f}%**\n')

    lines.append('\n## Summary Grid\n')
    if summary:
        rel = Path(summary).relative_to(RESOURCES.parent).as_posix()
        lines.append(f'![Summary]({rel})\n')

    lines.append('\n## Per-Model Detail\n')
    for case in cases:
        cd = ccx_all.get(case.name, {})
        ad = abq_all.get(case.name, {})
        lines.append(f'### {case.label}\n')
        c_ok = cd.get('success', False) if cd else False
        a_ok = ad.get('success', False) if ad else False
        lines.append(f'- CCX UMAT: {"OK" if c_ok else "FAIL"} (FLAG_HE={case.flag_he}, {len(case.ccx_umat_constants)} constants)')
        lines.append(f'- Abaqus: {"OK" if a_ok else "FAIL"} ({case.abq_he_keyword})')
        if cd and ad and cd.get('stress_11') and ad.get('stress_11'):
            cs = np.array(cd['stress_11'])
            a_s = np.array(ad['stress_11'])
            lines.append(f'- RMSE: {compute_rmse(cs, a_s):.3f}%')
            lines.append(f'- CCX last: {cs[-1]:.2f} MPa, ABQ last: {a_s[-1]:.2f} MPa')
        lines.append('')
        fname = RESOURCES / f'abq_ccx_{case.name}.png'
        if fname.exists():
            lines.append(f'![{case.name}]({fname.relative_to(RESOURCES.parent).as_posix()})\n')
        lines.append('')

    if avg < 2.0:
        lines.append('## Conclusion\n**Excellent agreement.** The CCX UMAT matches Abaqus built-in hyperelastic with average RMSE of '
                     f'**{avg:.3f}%**. The PRF UMAT implementation is verified against the industry-standard Abaqus solver.\n')
    elif avg < 5.0:
        lines.append('## Conclusion\n**Good agreement.** Average RMSE **{avg:.2f}%** — minor differences from volumetric formulation or tangent computation.\n')
    else:
        lines.append(f'## Conclusion\n**Partial agreement.** Average RMSE **{avg:.2f}%**. See per-model details above.\n')

    REPORT_PATH.write_text(''.join(lines), encoding='utf-8')
    return REPORT_PATH


# ─── Main ─────────────────────────────────────────────────────────────────

def main():
    print(f'Running Abaqus vs CalculiX comparison for {len(CASES)} models...\n')
    ccx_all: Dict[str, Dict] = {}
    abq_all: Dict[str, Dict] = {}

    for i, case in enumerate(CASES):
        print(f'[{i+1}/{len(CASES)}] {case.name}: {case.label}')

        ccx_inp = build_ccx_inp(case)
        abq_inp = build_abq_inp(case)

        print(f'  CCX...', end=' ')
        ccx_ok = run_ccx(ccx_inp.stem)
        print('[OK]' if ccx_ok else '[FAIL]')
        ccx_all[case.name] = parse_ccx(ccx_inp.stem) or {}

        print(f'  ABQ...', end=' ')
        abq_ok = run_abq(abq_inp.stem)
        print('[OK]' if abq_ok else '[FAIL]')
        abq_all[case.name] = parse_abq(abq_inp.stem) or {}

        print(f'  Plot...', end=' ')
        plot(case, ccx_all[case.name], abq_all[case.name])
        print('[OK]')

    print('\nSummary...', end=' ')
    sp = plot_summary(CASES, ccx_all, abq_all)
    print('[OK]')
    print('Report...', end=' ')
    r = generate_report(CASES, ccx_all, abq_all, str(sp))
    print('[OK]')
    print(f'\nReport: {r}')
    print('Done.')


if __name__ == '__main__':
    main()
