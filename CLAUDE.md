# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Does

AutoCalculix is a Python-based FEA (Finite Element Analysis) pipeline for modal (eigenvalue) analysis. It automates the full workflow: geometry → mesh → CalculiX solver → result parsing → 3D visualization. The primary use case is extracting natural frequencies and mode shapes from arbitrary Abaqus-compatible `.inp` files.

## Running the Pipeline

```bash
python src/pipeline.py
```

This runs the complete 6-step pipeline: mesh generation, model building, solving, DAT parsing, FRD→VTK conversion, and visualization.

## External Dependencies

- **CalculiX executable:** Hardcoded in `src/core/config.py` as `D:\SOFTWARE\calculix_2.23_4win\ccx_static.exe`
- **Python packages:** `gmsh`, `pyvista`, `numpy` (no requirements.txt — install manually)
- All output artifacts go to `workspace/`

## Architecture

`pipeline.py` is the orchestrator that calls six single-responsibility modules in sequence:

| Module | Class | Input → Output |
|--------|-------|----------------|
| `core/mesher.py` | `GmshMesher` | `TrayGeometryConfig` → `workspace/mesh.inp` |
| `core/model_builder.py` | `CalculixModelBuilder` | `mesh.inp` + `ModalAnalysisConfig` → `workspace/<job>.inp` (master INP that `*INCLUDE`s mesh) |
| `core/solver.py` | `CalculixSolver` | `<job>.inp` → `<job>.frd`, `<job>.dat` |
| `core/dat_parser.py` | `CalculixDatParser` | `<job>.dat` → eigenfrequencies list (Hz) |
| `core/frd_converter.py` | `FrdToVtuConverter` | `<job>.frd` → `<job>.01.vtu`, `<job>.02.vtu`, … (ccx2paraview, 1-based 2-digit) |
| `core/viewer.py` | `ModeShapeViewer` | `<job>.{N:02d}.vtu` → warped PyVista window |

**Configuration** for paths, material properties (Steel: E=210GPa, ν=0.3, ρ=7.85e-9 t/mm³), and modal settings (default 10 modes) lives in `src/core/config.py`.

## Key Technical Details

- **INP structure:** Master INP uses `*INCLUDE` to reference the mesh INP — this pattern supports multi-part geometries.
- **Element conversion:** Gmsh outputs 2D plane elements (CPS/CPE); `mesher.py` converts these to CalculiX shell elements (S3/S4).
- **FRD → VTU:** Uses `ccx2paraview` (`pip install ccx2paraview`). Outputs one VTU per mode: `<job>.01.vtu`, `<job>.02.vtu`, … (1-based, 2-digit zero-padded).
- **Shell element expansion:** CalculiX internally upgrades S3 → wedge elements with z=±thickness/2 nodes in FRD output.
- **Rigid body modes:** Free-free (unconstrained) analysis yields ~7 rigid body / near-zero modes; the first meaningful flexible mode is typically mode 8.
- **Displacement scaling in viewer:** Warped geometry is scaled to 10% of the maximum model dimension for visibility.
- **Two pipeline entry points:** `run_with_meshing(analysis, geometry)` for Gmsh-generated mesh; `run_from_inp(mesh_inp_file, analysis)` for existing INP files.
