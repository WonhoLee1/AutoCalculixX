from pathlib import Path
p = Path('D:\\PythonCodeStudy\\AutoCalculix\\workspace\\plot_n0.frd')
txt = p.read_text(encoding='utf-8', errors='replace')
lines = txt.splitlines()
print(f"Total lines: {len(lines)}")

in_stress = False
in_strain = False
stress_count = 0
strain_count = 0

for i, line in enumerate(lines):
    parts = line.strip().split()
    if not parts:
        continue
    if len(parts) >= 3 and parts[0] == '-4' and parts[1] == 'STRESS':
        in_stress = True
        in_strain = False
        print(f"L{i}: STRESS block start")
        continue
    if len(parts) >= 3 and parts[0] == '-4' and parts[1] == 'TOSTRAIN':
        in_strain = True
        in_stress = False
        print(f"L{i}: TOSTRAIN block start")
        continue
    if parts[0] == '-3':
        if in_stress: print(f"L{i}: STRESS end (was {stress_count} pts)")
        in_stress = False
        in_strain = False
        continue
    if parts[0] == '-5':
        continue
    if parts[0] == '-1':
        if in_stress:
            stress_count += 1
            if stress_count == 1:
                print(f"L{i}: FIRST stress -1 line: parts={parts}, vals={parts[2:8]}")
        if in_strain:
            strain_count += 1
            if strain_count == 1:
                print(f"L{i}: FIRST strain -1 line: parts={parts}, vals={parts[2:8]}")
        continue

print(f"\nTotal stress pts: {stress_count}, strain pts: {strain_count}")
