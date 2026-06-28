# WHTOOLs MATCALIB 2026 — Reproducibility test
import sys, subprocess, os
sys.path.insert(0, 'D:\\PythonCodeStudy\\AutoCalculix')
from src.core.config import PRFConfig, ViscousNetworkConfig, PolynomialHyperelasticParams
from src.core.model_builder import CalculixModelBuilder
from pathlib import Path

os.environ['Path'] = r'C:\msys64\mingw64\bin;' + os.environ.get('Path', '')
CCX = r'D:\SOFTWARE\calculix_2.23_4win\ccx_custom.exe'
WS = Path('D:\\PythonCodeStudy\\AutoCalculix\\workspace')

print('Building INP...', flush=True)
cfg = PRFConfig(
    job_name='test_repro', flag_he=1,
    he_poly=PolynomialHyperelasticParams(C10=10.0, C01=0.0, RBULK_HE=100.0),
    networks=[ViscousNetworkConfig(stiffn=0.5, A1=0.1, EXPC=-1.0, EXPM=2.0, KSI=0.01)],
    initial_inc=0.01, step_time=1.0, max_inc=0.1,
)
# clean first
for f in WS.glob('test_repro.*'):
    f.unlink(missing_ok=True)

builder = CalculixModelBuilder(WS)
inp = builder.build_prf_inp(cfg)
print(f'INP: {inp}', flush=True)

print('Running CCX...', flush=True)
result = subprocess.run([CCX, inp.stem], cwd=str(WS), capture_output=True, text=True, timeout=300)
print(f'returncode: {result.returncode}', flush=True)
print('stdout:', result.stdout[-500:], flush=True)
print('stderr:', result.stderr[-500:], flush=True)

sta = WS / 'test_repro.sta'
if sta.exists():
    lines = [l.strip() for l in sta.read_text().splitlines() if l.strip() and not l.startswith('SUMMARY') and 'STEP' not in l and 'INC' not in l and 'ATT' not in l]
    print(f'Total incs: {len(lines)}', flush=True)
    if lines:
        print(f'First: {lines[0]}', flush=True)
        print(f'Last: {lines[-1]}', flush=True)
else:
    print('STA file not found', flush=True)

print('DONE', flush=True)
