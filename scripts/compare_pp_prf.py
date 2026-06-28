#!/usr/bin/env python3
"""PP C3TF2 PRF Model: CCX UMAT vs Experimental Data"""
import sys, os, re, subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ['Path'] = r'C:\msys64\mingw64\bin;' + os.environ.get('Path', '')

CCX_EXE = r'D:\SOFTWARE\calculix_2.23_4win\ccx_custom.exe'
DATA_DIR = ROOT.parent/'WHT_PRF'/'doc'/'Calibration of a PRF Material Model for Polypropylene'/'PPC3TF2TestData'
WORKSPACE = ROOT/'workspace'
RESOURCES = ROOT/'dev_log'/'resources'
REPORT_PATH = ROOT/'dev_log'/'pp_prf_validation.md'
RESOURCES.mkdir(parents=True, exist_ok=True)
WORKSPACE.mkdir(parents=True, exist_ok=True)

# ── Material Model ───────────────────────────────────────────────────────────
C10 = 549.6
mu = 2 * C10
nu = 0.45
K = 2 * mu * (1 + nu) / (3 * (1 - 2 * nu))
D1_CCX = K / 2.0
print(f"C10={C10}, mu={mu:.0f}, nu={nu}, K={K:.0f}, D1={D1_CCX:.0f}")

CCX_CONSTANTS = [
    3.0, 1.0, 0.0,           # N_NET=3, FLAG_HE=1, FLAG_PL=0
    C10, 0,0,0,0, 0,0,0,0,   # C10..C03
    D1_CCX, 0.0, 0.0,         # D1, D2, D3
    1.0, K,                   # IHYPER=1, RBULK
    0.337060, 1.0, 4.8828e-7, -0.551660, 3.13398, 0.001, 1.0,  # NW1: BB
    0.150654, 1.0, 2.4414e-6, -0.746191, 4.91992, 0.001, 1.0,  # NW2: BB
    0.372607, 1.0, 7.81403e-5, -0.616699, 3.801953, 0.001, 1.0, # NW3: BB
    1.0, K,                   # G, RBULK_FINAL
]
N_CONSTANTS = len(CCX_CONSTANTS)
DEPVAR = 2 + 12 * 3
print(f"Constants: {N_CONSTANTS}, DEPVAR: {DEPVAR}")

MESH = """*NODE
1,0,0,0\n2,1,0,0\n3,1,1,0\n4,0,1,0\n5,0,0,1\n6,1,0,1\n7,1,1,1\n8,0,1,1
*ELEMENT,TYPE=C3D8,ELSET=EA\n1,1,2,3,4,5,6,7,8\n"""

BC_FIX = "1,1,3\n4,1,1\n4,3,3\n8,1,1\n5,1,2\n2,2,3\n6,2,2\n2,3,3\n3,3,3"
BC_PULL = lambda d: f"2,1,1,{d:.10f}\n3,1,1,{d:.10f}\n6,1,1,{d:.10f}\n7,1,1,{d:.10f}"

_FLOAT_RE = re.compile(r'[-+]?\d+\.\d+E[+-]\d+')

# ── INP builder ─────────────────────────────────────────────────────────────

def fmt(c): 
    a = list(c) + [0.0]
    return ''.join(', '.join(f'{v:.15g}' for v in a[i:i+8]) + '\n' for i in range(0,len(a),8) if a[i:i+8])

def build_inp(name, disp, t, inc, material_section=None):
    out = WORKSPACE / f'{name}.inp'
    mat = material_section or f"*SOLID SECTION,ELSET=EA,MATERIAL=USERPRF\n1.\n\n*MATERIAL,NAME=USERPRF\n*USER MATERIAL,CONSTANTS={N_CONSTANTS}\n{fmt(CCX_CONSTANTS)}\n*DEPVAR\n{DEPVAR}\n"
    out.write_text(f"*HEADING\n{name}\n{MESH}\n{mat}\n*STEP,INC=5000,NLGEOM=YES\n*STATIC\n{inc},{t},1e-20,{t/10}\n*BOUNDARY\n{BC_FIX}\n{BC_PULL(disp)}\n*EL FILE\nS,E\n*END STEP\n")
    return out

def build_relax_inp(name, max_disp, ramp_t, hold_t):
    out = WORKSPACE / f'{name}.inp'
    mat = f"*SOLID SECTION,ELSET=EA,MATERIAL=USERPRF\n1.\n\n*MATERIAL,NAME=USERPRF\n*USER MATERIAL,CONSTANTS={N_CONSTANTS}\n{fmt(CCX_CONSTANTS)}\n*DEPVAR\n{DEPVAR}\n"
    parts = [f"*HEADING\n{name}\n{MESH}\n{mat}"]
    for step_t, label in [(ramp_t, "RAMP"), (hold_t, "HOLD")]:
        parts.append(f"*STEP,INC=500,NLGEOM=YES\n*STATIC\n0.01,{step_t},1e-20,{step_t/10}\n*BOUNDARY\n{BC_FIX}\n{BC_PULL(max_disp)}\n*EL FILE\nS,E\n*END STEP\n")
    out.write_text(''.join(parts))
    return out

# ── Runner & Parser ─────────────────────────────────────────────────────────

def run_ccx(stem): 
    try: subprocess.run([CCX_EXE, stem], cwd=str(WORKSPACE), capture_output=True, text=True, timeout=600); return (WORKSPACE/f'{stem}.frd').exists()
    except: return False

def parse_frd(stem):
    frd = WORKSPACE / f'{stem}.frd'
    if not frd.exists(): return {'stress_11':[], 'strain_11':[], 'success':False}
    text = frd.read_text(encoding='utf-8', errors='replace')
    r={'stress_11':[],'strain_11':[],'success':False}; sm=em=dm=False
    for line in text.splitlines():
        s=line.strip()
        if not s: continue
        p=s.split()
        if p[0]=='-4':
            if len(p)>=2 and p[1]=='STRESS': sm,em,dm=True,False,False
            elif len(p)>=2 and p[1]=='TOSTRAIN': em,sm,dm=True,False,False
            continue
        if p[0] in('-3','-5'): continue
        if p[0]=='-1' and (sm or em) and not dm and len(p)>=2 and p[1]=='1':
            fl=[float(m) for m in _FLOAT_RE.findall(line)]
            if len(fl)>=2: (r['stress_11'] if sm else r['strain_11']).append(fl[0]); dm=True
    sta = WORKSPACE/f'{stem}.sta'
    if sta.exists():
        for line in sta.read_text(encoding='utf-8',errors='replace').splitlines():
            s=line.strip()
            if not s or 'SUMMARY' in s: continue
            parts=re.split(r'\s+',s)
            if len(parts)>=5:
                try:
                    if float(parts[4].replace('U',''))>0.95: r['success']=True
                except: pass
    return r

# ── Load experimental data ──────────────────────────────────────────────────

rate_exp = np.loadtxt(DATA_DIR/'rate_data_1e+2_partial.txt')
relax_exp = {p: np.loadtxt(DATA_DIR/f'relax_data_{p}.txt', delimiter=',', comments='#') for p in ['0050','0075','0100','0150']}
print(f"Rate: {rate_exp.shape}, Relax: " + ", ".join(f"{k}:{v.shape[0]}" for k,v in relax_exp.items()))

# ── Rate Test ───────────────────────────────────────────────────────────────

print("\n=== Rate Test ===")
inp = build_inp('pp_rate', rate_exp[-1,1], 1.0, 0.02)
print(f"  INP: {inp.name}, max_disp={rate_exp[-1,1]:.6f}")
ok=run_ccx(inp.stem); print(f"  CCX: {'OK' if ok else 'FAIL'}")
rate_ccx=parse_frd(inp.stem)
print(f"  Stress pts: {len(rate_ccx['stress_11'])}, Success: {rate_ccx['success']}")

# ── Relaxation Tests ────────────────────────────────────────────────────────

relax_ccx = {}
for pct, edata in relax_exp.items():
    print(f"\n=== Relax {pct} ===")
    max_strain = edata[-1,1]
    hold_time = edata[-1,0]
    ramp_time = hold_time * 0.01
    inp = build_relax_inp(f'pp_relax_{pct}', max_strain, ramp_time, hold_time)
    print(f"  INP: {inp.name}, strain={max_strain:.6f}, ramp={ramp_time:.1f}s, hold={hold_time:.0f}s")
    ok=run_ccx(inp.stem); print(f"  CCX: {'OK' if ok else 'FAIL'}")
    relax_ccx[pct]=parse_frd(inp.stem)
    print(f"  Stress pts: {len(relax_ccx[pct]['stress_11'])}, Success: {relax_ccx[pct]['success']}")

# ── Plotting ────────────────────────────────────────────────────────────────

EC,CC='#333','#2196F3'

# Rate test
fig,ax=plt.subplots(figsize=(9,5))
ax.plot(rate_exp[:,1]*100, rate_exp[:,2], color=EC, lw=2, label='Experiment')
if rate_ccx['stress_11']:
    s=np.array(rate_ccx['stress_11']); e=np.array(rate_ccx['strain_11'][:len(s)])
    ax.plot(e*100, s, color=CC, lw=2, ls='--', marker='.', ms=2, label='CCX UMAT')
ax.set_xlabel('Engineering Strain (%)'); ax.set_ylabel('Engineering Stress (MPa)')
ax.set_title('PP C3TF2 — Rate Test 100/s'); ax.legend(); ax.grid(True,alpha=.3)
f=RESOURCES/'pp_rate_100.png'; fig.savefig(f,dpi=150,bbox_inches='tight'); plt.close(fig)

# Relaxation grid
fig,axes=plt.subplots(2,2,figsize=(12,9)); axes=axes.flatten()
for idx,(pct,edata) in enumerate(relax_exp.items()):
    ax=axes[idx]
    ax.plot(edata[:,0], edata[:,2], color=EC, lw=1.5, label=f'Exp {pct}%')
    rd=relax_ccx.get(pct)
    if rd and rd['stress_11']:
        s=np.array(rd['stress_11']); t=np.linspace(0,edata[-1,0],len(s))
        ax.plot(t, s, color=CC, lw=2, ls='--', label='CCX UMAT')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Eng. Stress (MPa)')
    ax.set_title(f'Relaxation {pct}% Strain'); ax.legend(fontsize=8); ax.grid(True,alpha=.3)
fig.suptitle('PP C3TF2 — Stress Relaxation: Experiment vs CCX UMAT',fontsize=13,fontweight='bold')
fig.tight_layout()
f=RESOURCES/'pp_relax_grid.png'; fig.savefig(f,dpi=150,bbox_inches='tight'); plt.close(fig)

# ── Report ──────────────────────────────────────────────────────────────────

import datetime
ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
REPORT_PATH.write_text(f"""# WHTOOLs MATCALIB 2026 — PP C3TF2 PRF Model Validation

**Generated**: {ts}

## Material Model

Abaqus calibrated PRF model (polypropylene PP C3TF2):

| Parameter | Value |
|---|---|
| Hyperelastic | Neo-Hookean, C10=549.6 MPa |
| Bulk modulus | K={K:.0f} MPa (ν={nu}) |
| NW1 (SRatio=0.337) | Law=strain, A=4.88e-7, m=3.13, C=-0.55 |
| NW2 (SRatio=0.151) | Law=strain, A=2.44e-6, m=4.92, C=-0.75 |
| NW3 (SRatio=0.373) | Law=strain, A=7.81e-5, m=3.80, C=-0.62 |
| CCX UMAT | N_NETWORK=3, FLAG_HE=1, {N_CONSTANTS} constants, DEPVAR={DEPVAR} |

## Rate Test (100/s)

![Rate Test](resources/pp_rate_100.png)

## Relaxation Tests

![Relaxation Grid](resources/pp_relax_grid.png)
""")
print(f"\nReport: {REPORT_PATH}")
print("Done.")
