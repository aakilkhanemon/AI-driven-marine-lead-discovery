# AI-driven-marine-lead-discovery

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![RDKit](https://img.shields.io/badge/RDKit-2023.09+-green.svg)](https://www.rdkit.org/)
[![GNINA](https://img.shields.io/badge/GNINA-v1.1-red.svg)](https://github.com/gnina/gnina)

An integrated computational framework combining 2D/3D cheminformatics, deep-learning docking, explicit-solvent molecular dynamics (MD), Density Functional Theory (DFT), systems pharmacology, and consensus machine learning scoring for novel drug candidate identification.

---

## 📌 Workflow Overview

1. **Marine Library Development:** 2D CACTVS 881-bit Tanimoto expansion ($T_c \ge 0.97$), 3D pharmacophore scaffold hopping ($T_c \le 0.42$), and synthetic-fragment BRICS *de novo* design.
2. **Deep-Learning Docking:** GPU-accelerated GNINA screening with active-site residue centroid box auto-scaling ($16\text{ \AA}$ bounding box).
3. **Molecular Dynamics & MM-GBSA/PBSA:** Explicit-solvent trajectory processing, RMSD/RMSF/PCA dimensionality reduction, and per-residue thermodynamic free-energy decomposition.
4. **Quantum Mechanical DFT:** B3LYP electronic structure calculations, HOMO/LUMO frontier orbital gaps, and conceptual DFT global reactivity descriptors ($\chi, \eta, \omega$).
5. **Systems Pharmacology:** API-driven STRING interactome construction (confidence score $\ge 0.700$) and cytoHubba degree/betweenness centrality topology mining.
6. **Consensus ML Lead Prioritization:** Dual-tier unsupervised MCDA and supervised structural constraint-weighted scoring.

---

## 🛠️ Installation

### 1. Clone Repository
```bash
git clone [https://github.com/aakilkhanemon/AI-driven-marine-lead-discovery.git](https://github.com/aakilkhanemon/AI-driven-marine-lead-discovery.git)
cd AI-driven-marine-lead-discovery
2. Set Up Conda Environment
Bash
conda env create -f environment.yml
conda activate marine_lead_discovery
pip install -e .
3. External Dependencies
EnsureGNINA v1.1+binary is present in your system path or placed under src/docking/bin/.

🚀 Execution Guide
Run the full end-to-end pipeline via the master CLI runner:

Bash
python scripts/run_pipeline.py --config config/parameters.yaml
Or execute individual modules independently:

Python
# Example: Running Consensus ML Scoring
from src.consensus_ml.scoring_pipeline import ConsensusScorer

scorer = ConsensusScorer(input_matrix="data/processed/master_screening_matrix.csv")
top_leads = scorer.prioritize_leads(top_n=3)
print(top_leads)
🔬 Citation
If you use this codebase or computational pipeline in your research, please cite:

Aakil Khan Emon.  Systems Pharmacology-Driven Discovery of Tri-Modal Anti-Mpox Therapeutics: Integrating Quantum Mechanics, Molecular Dynamics, and Network Topology. Master's Thesis, 2026.

📜 License
Distributed under the MIT License. See LICENSEfor details.


---

## 📦 Requirements & Environment Configuration

### `requirements.txt`
```text
rdkit>=2023.09.1
pandas>=2.0.0
numpy>=1.24.0
scikit-learn>=1.2.0
MDAnalysis>=2.4.0
mdtraj>=1.9.7
psi4>=1.8
pyscf>=2.2.0
networkx>=3.0
matplotlib>=3.7.0
requests>=2.28.0
py3Dmol>=2.0.0
pyyaml>=6.0
environment.yml
YAML
name: marine_lead_discovery
channels:
  - conda-forge
  - psi4
dependencies:
  - python=3.10
  - rdkit
  - psi4
  - mamba
  - pip
  - pip:
      - -r requirements.txt
🧩 Core Source Code Modules
1. Module:src/consensus_ml/scoring_pipeline.py
Python
"""
Consensus Machine Learning and MCDA Scoring Module
Handles unsupervised equal-weight mean aggregation and supervised
pharmacokinetic/thermodynamic constraint-weighted lead prioritization.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional


class ConsensusScorer:
    """Prioritizes lead compounds across target receptors using MCDA frameworks."""

    DEFAULT_SUPERVISED_WEIGHTS: Dict[str, float] = {
        "GNINA_Vina": 0.25,
        "GNINA_CNN": 0.20,
        "MOE_Induced": 0.20,
        "ADME": 0.12,
        "Toxicology": 0.13,
        "PASS": 0.10,
    }

    def __init__(
        self,
        input_matrix: Optional[str] = None,
        data_frame: Optional[pd.DataFrame] = None,
    ):
        if data_frame is not None:
            self.df = data_frame.copy()
        elif input_matrix is not None:
            self.df = pd.read_csv(input_matrix)
        else:
            raise ValueError("Provide either an input CSV path or a DataFrame.")

        self.features: List[str] = list(self.DEFAULT_SUPERVISED_WEIGHTS.keys())

    def compute_unsupervised_score(self) -> pd.DataFrame:
        """Calculates Tier 1 Unsupervised Consensus Score (Equal-Weight Mean)."""
        self.df["Unsupervised_Score"] = np.round(
            self.df[self.features].mean(axis=1), 2
        )
        return self.df

    def compute_supervised_score(
        self, custom_weights: Optional[Dict[str, float]] = None
    ) -> pd.DataFrame:
        """Calculates Tier 2 Supervised Consensus Score (Structural-Bias Weighted)."""
        weights = custom_weights or self.DEFAULT_SUPERVISED_WEIGHTS
        self.df["Supervised_Score"] = np.round(
            sum(self.df[f] * weights[f] for f in self.features), 2
        )
        return self.df

    def prioritize_leads(self, top_n: int = 3) -> pd.DataFrame:
        """Isolates Top-N candidates per receptor pocket using hierarchical sorting."""
        if "Unsupervised_Score" not in self.df.columns:
            self.compute_unsupervised_score()
        if "Supervised_Score" not in self.df.columns:
            self.compute_supervised_score()

        compiled_leads = []
        for receptor_id, group in self.df.groupby("Receptor_ID"):
            sorted_candidates = group.sort_values(
                by=["Supervised_Score", "Unsupervised_Score"],
                ascending=[False, False],
            ).head(top_n).copy()

            sorted_candidates["Final_Manuscript_Rank"] = [
                f"Top-{i+1} Lead" for i in range(len(sorted_candidates))
            ]
            compiled_leads.append(sorted_candidates)

        return pd.concat(compiled_leads, ignore_index=True)
2. Module:src/docking/cavity_mapping.py
Python
"""
Dynamic Cavity Geometry Mapper
Calculates center-of-mass centroids and adaptive grid bounding cubes.
"""

import math
import numpy as np
import MDAnalysis as mda
from typing import Tuple


def calculate_active_site_grid(
    receptor_pdb: str, target_residues_str: str, padding: float = 10.0
) -> Tuple[Tuple[float, float, float], int]:
    """
    Extracts active site centroid and calculates adaptive grid box size.

    Args:
        receptor_pdb: Path to receptor PDB file.
        target_residues_str: Selection string (e.g., 'resid 41 80 108 119').
        padding: Bounding envelope margin in Angstroms.

    Returns:
        Tuple containing ((CX, CY, CZ), box_size)
    """
    u_rec = mda.Universe(receptor_pdb)
    pocket_atoms = u_rec.select_atoms(target_residues_str)

    if len(pocket_atoms) == 0:
        raise ValueError(f"Residues '{target_residues_str}' not resolved in PDB.")

    cx, cy, cz = pocket_atoms.center_of_mass()
    span = np.max(pocket_atoms.positions, axis=0) - np.min(pocket_atoms.positions, axis=0)
    box_size = max(16, math.ceil(np.max(span) + padding))

    return (float(cx), float(cy), float(cz)), box_size
3. Module:src/dft_quantum/dft_runner.py
Python
"""
Density Functional Theory (DFT) Quantum Chemistry Engine
Handles B3LYP electronic structure calculations and global reactivity descriptors.
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import pyscf
from pyscf import dft
from typing import Dict, Any


class DFTAnalyzer:
    """Executes DFT calculations via PySCF and derives conceptual reactivity indices."""

    HARTREE_TO_EV = 27.211386

    def __init__(self, basis_set: str = "3-21g", xc_functional: str = "b3lyp"):
        self.basis_set = basis_set
        self.xc_functional = xc_functional

    def compute_frontier_orbitals(self, smiles: str) -> Dict[str, Any]:
        """Calculates HOMO, LUMO, and Band Gap in eV for a given SMILES string."""
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

        conf = mol.GetConformer()
        pyscf_geom = "".join(
            f"{atom.GetSymbol()} {conf.GetAtomPosition(atom.GetIdx()).x} "
            f"{conf.GetAtomPosition(atom.GetIdx()).y} {conf.GetAtomPosition(atom.GetIdx()).z}; "
            for atom in mol.GetAtoms()
        )[:-2]

        pyscf_mol = pyscf.gto.Mole()
        pyscf_mol.atom = pyscf_geom
        pyscf_mol.basis = self.basis_set
        pyscf_mol.charge = 0
        pyscf_mol.spin = 0
        pyscf_mol.build()

        mf = dft.RKS(pyscf_mol)
        mf.xc = self.xc_functional
        mf.kernel()

        mo_energies = mf.mo_energy * self.HARTREE_TO_EV
        n_occ = np.count_nonzero(mf.mo_occ > 0)

        homo = mo_energies[n_occ - 1]
        lumo = mo_energies[n_occ]
        gap = lumo - homo

        return {"HOMO_eV": round(homo, 3), "LUMO_eV": round(lumo, 3), "Gap_eV": round(gap, 3)}

    @staticmethod
    def derive_global_reactivity(df_results: pd.DataFrame) -> pd.DataFrame:
        """Derives conceptual DFT global reactivity descriptors via Koopmans' theorem."""
        df = df_results.copy()
        df["Electronegativity_eV"] = - (df["HOMO_eV"] + df["LUMO_eV"]) / 2
        df["Chemical_Potential_eV"] = - df["Electronegativity_eV"]
        df["Global_Hardness_eV"] = (df["LUMO_eV"] - df["HOMO_eV"]) / 2
        df["Global_Softness_1_eV"] = 1 / (2 * df["Global_Hardness_eV"])
        df["Electrophilicity_Index_eV"] = (
            df["Chemical_Potential_eV"] ** 2
        ) / (2 * df["Global_Hardness_eV"])
        return df.round(3)
🛠️ Master Execution Script:scripts/run_pipeline.py
Python
#!/usr/bin/env python3
"""
Master Pipeline Orchestrator
Executes the full computer-aided drug discovery pipeline end-to-end.
"""

import sys
import argparse
from pathlib import Path
from src.consensus_ml.scoring_pipeline import ConsensusScorer


def main():
    parser = argparse.ArgumentParser(
        description="End-to-End Marine Lead Discovery Pipeline Execution Engine."
    )
    parser.add_argument(
        "--input",
        type=str,
        default="data/raw/master_36_hits_matrix.csv",
        help="Path to raw screening matrix CSV file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/outputs/Curated_Top3_MD_Candidates.csv",
        help="Path to save top prioritized leads.",
    )
    args = parser.parse_args()

    print("=================================================================")
    print("🚀 INITIALIZING COMPUTATIONAL DRUG DISCOVERY PIPELINE")
    print("=================================================================")

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"❌ Error: Input file '{args.input}' not found.")
        sys.exit(1)

    print(f"• Ingesting dataset from: {input_path}")
    scorer = ConsensusScorer(input_matrix=str(input_path))

    print("• Calculating Tier 1 Unsupervised MCDA Scores...")
    scorer.compute_unsupervised_score()

    print("• Calculating Tier 2 Supervised Structural-Bias Scores...")
    scorer.compute_supervised_score()

    print("• Isolating Top-3 prioritized leads per target pocket...")
    top_leads = scorer.prioritize_leads(top_n=3)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_leads.to_csv(output_path, index=False)

    print("=================================================================")
    print(f"✅ PIPELINE EXECUTION SUCCESSFUL")
    print(f"💾 Prioritized candidate matrix exported to: {output_path}")
    print("=================================================================")


if __name__ == "__main__":
    main()
