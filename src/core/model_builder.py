# -*- coding: utf-8 -*-
# WHTOOLs MATCALIB 2026 — CalculiX model builder
from pathlib import Path
from src.core.config import (
    ModalAnalysisConfig, HyperelasticConfig, PRFConfig,
    VISC_BB, VISC_SINH, VISC_POWER,
)


class CalculixModelBuilder:
    def __init__(self, workspace: Path):
        self.workspace = workspace
        self.workspace.mkdir(parents=True, exist_ok=True)

    def build_modal_analysis(self, mesh_inp_file: Path, config: ModalAnalysisConfig) -> Path:
        """
        메쉬 INP를 *INCLUDE하는 Master INP를 생성합니다.
        mesh_inp_file은 workspace 안에 있어야 합니다 (CalculiX가 상대경로로 INCLUDE).
        """
        master_inp_path = self.workspace / f"{config.job_name}.inp"

        with open(master_inp_path, 'w', encoding='utf-8') as f:
            f.write(f"*HEADING\nModal Analysis: {config.job_name}\n")
            f.write(f"*INCLUDE, INPUT={mesh_inp_file.name}\n\n")

            f.write(f"*MATERIAL, NAME=MAT1\n")
            f.write(f"*ELASTIC\n{config.E}, {config.nu}\n")
            f.write(f"*DENSITY\n{config.rho}\n\n")

            f.write(f"*SHELL SECTION, ELSET={config.elset_name}, MATERIAL=MAT1\n")
            f.write(f"{config.thickness}\n\n")

            f.write("*STEP\n*FREQUENCY\n")
            f.write(f"{config.num_modes}\n")
            f.write("*NODE FILE\nU\n")
            f.write("*END STEP\n")

        return master_inp_path

    def build_hyperelastic_inp(self, config: HyperelasticConfig) -> Path:
        """
        Creates a single-element *USER MATERIAL INP file for
        polynomial hyperelasticity with NLGEOM static step.

        The *USER MATERIAL card writes 14 constants in the exact
        order expected by src/umat/hyperelastic_user.f.
        """
        inp = self.workspace / f"{config.job_name}.inp"

        coords = [
            "0.0, 0.0, 0.0",
            "1.0, 0.0, 0.0",
            "1.0, 1.0, 0.0",
            "0.0, 1.0, 0.0",
            "0.0, 0.0, 1.0",
            "1.0, 0.0, 1.0",
            "1.0, 1.0, 1.0",
            "0.0, 1.0, 1.0",
        ]

        lines = [
            f"*HEADING\nPolynomial Hyperelastic: {config.job_name}\n",
            "*NODE\n",
        ]
        for i, xyz in enumerate(coords, 1):
            lines.append(f"{i}, {xyz}\n")

        nids = list(range(1, 9))
        lines.append("*ELEMENT,TYPE=C3D8,ELSET=EL1\n")
        lines.append(f"1, {', '.join(map(str, nids))}\n\n")
        lines.append(f"*SOLID SECTION,ELSET=EL1,MATERIAL={config.material_name}\n\n")
        lines.append(f"*MATERIAL,NAME={config.material_name}\n")
        lines.append("*USER MATERIAL,CONSTANTS=14\n")

        c = list(config.constants)
        c.append(0.0)  # temperature (elcon(0))
        for i in range(0, 15, 8):
            chunk = c[i:i+8]
            while len(chunk) < 8:
                chunk.append(0.0)
            lines.append(", ".join(f"{v:.15g}" for v in chunk) + "\n")
        lines.append("\n")
        lines.append("*STEP,INC=100,NLGEOM=YES\n")
        lines.append("*STATIC\n")
        lines.append(f"{config.initial_inc}, {config.step_time}, "
                     f"{config.min_inc}, {config.max_inc}\n")
        lines.append("*BOUNDARY\n")

        bcs = [
            ("1,1,3\n"),
            ("2,2,3\n"),
            ("4,2,3\n"),
            ("5,1,3\n"),
            ("6,2,3\n"),
            ("8,2,3\n"),
            ("1,1,1\n"),
            ("2,1,1,0.1\n"),
            ("3,1,1,0.1\n"),
            ("4,1,1\n"),
            ("5,1,1\n"),
            ("6,1,1,0.1\n"),
            ("7,1,1,0.1\n"),
            ("8,1,1\n"),
        ]
        lines.extend(bcs)
        lines.append("*NODE FILE\nU\n")
        lines.append("*EL FILE\nS,E\n")
        lines.append("*END STEP\n")

        with open(inp, 'w', encoding='utf-8') as f:
            f.writelines(lines)

        print(f"[Builder] Created hyperelastic INP: {inp}")
        return inp

    def build_prf_inp(self, config: PRFConfig) -> Path:
        """
        Creates a single-element PRF (Multi-Network Foam) UMAT INP.
        Supports all viscous models: BB, Sinh, Power.
        Uses C3D8 with NLGEOM static step + creep.
        """
        inp = self.workspace / f"{config.job_name}.inp"

        n_nodes = 8
        coords = [
            "0.0, 0.0, 0.0",
            "1.0, 0.0, 0.0",
            "1.0, 1.0, 0.0",
            "0.0, 1.0, 0.0",
            "0.0, 0.0, 1.0",
            "1.0, 0.0, 1.0",
            "1.0, 1.0, 1.0",
            "0.0, 1.0, 1.0",
        ]

        lines = [
            f"*HEADING\nPRF Multinetwork: {config.job_name}\n",
            "*NODE\n",
        ]
        for i, xyz in enumerate(coords, 1):
            lines.append(f"{i}, {xyz}\n")

        nids = list(range(1, n_nodes + 1))
        lines.append(f"*ELEMENT,TYPE={config.element_type},ELSET=EL1\n")
        lines.append(f"1, {', '.join(map(str, nids))}\n\n")
        lines.append(f"*SOLID SECTION,ELSET=EL1,MATERIAL={config.material_name}\n\n")
        lines.append(f"*MATERIAL,NAME={config.material_name}\n")

        nc = config.nconstants
        lines.append(f"*USER MATERIAL,CONSTANTS={nc}\n")
        cons = list(config.constants)
        cons.append(0.0)  # elcon(0) = temperature reference
        for i in range(0, nc + 1, 8):
            chunk = cons[i:i+8]
            if not chunk:
                continue
            # Pad last line with zeros to exactly 8 values (CalculiX parser)
            while len(chunk) < 8:
                chunk.append(0.0)
            lines.append(", ".join(f"{v:.15g}" for v in chunk) + "\n")
        lines.append("\n")

        lines.append(f"*DEPVAR\n{config.depvar}\n\n")
        lines.append("*STEP,INC=200,NLGEOM=YES\n")
        lines.append("*STATIC\n")
        lines.append(f"{config.initial_inc}, {config.step_time}, "
                     f"{config.min_inc}, {config.max_inc}\n")
        lines.append("*BOUNDARY\n")

        bcs = [
            "1,1,3\n",
            "2,2,3\n",
            "4,2,3\n",
            "5,1,3\n",
            "6,2,3\n",
            "8,2,3\n",
            "1,1,1\n",
            "2,1,1,0.1\n",
            "3,1,1,0.1\n",
            "4,1,1\n",
            "5,1,1\n",
            "6,1,1,0.1\n",
            "7,1,1,0.1\n",
            "8,1,1\n",
        ]
        lines.extend(bcs)
        lines.append("*NODE FILE\nU\n")
        lines.append("*EL FILE\nS,E\n")
        lines.append("*END STEP\n")

        with open(inp, 'w', encoding='utf-8') as f:
            f.writelines(lines)

        print(f"[Builder] Created PRF INP: {inp}  ({nc} constants, {config.depvar} depvars)")
        return inp
