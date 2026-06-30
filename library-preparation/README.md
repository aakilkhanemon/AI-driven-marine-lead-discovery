The Defined Pre-Docking Pipeline (Chronological Order)
Phase I: The 2D Filter Gates (Max Library Reduction)
1.	Deduplication: Canonicalize SMILES and remove redundant structures from the initial massive-compound pool to establish a unique baseline dataset.
2.	Structure Sanitization: Evaluate basic chemical valency, standardize aromaticity models, and rectify structural inconsistencies.
3.	Fragment Removal / Desalting: Strip counter-ions, solvents, and salt placeholders (Na+, Cl-, etc.) to isolate the core active marine scaffold.
4.	Functional Group Standardization: Normalize volatile or variable representations of complex functional arrangements (eg, nitro, azide, or sulfonyl groups).
5.	Chemoinformatic discard Bioavailability Filtering: Apply hard property cuts (Lipinski's Rule of 5, Veber's criteria, REOS) while the library is still flat to non-drug-like space.
6.	Pan-Assay Interference Compounds (PAINS) Filtering: Screen out highly reactive, aggregative, or assay-disrupting substructures prior to any 3D matrix scaling.
Phase II: Physicochemical and Stereochemical Expansion
7.	Physiological Protonation & Ionization (pH 7.4): Calculate dominant protonation, tautomerization, and macrocyclic ionization states strictly under simulated biological conditions.
8.	Stereochemical Fidelity & Enumeration: Generate explicit 3D representations for unassigned chiral centers, undefined double bonds, or flexible macrocyclic cis/trans geometry.
Phase III: 3D Generation and Quantum Parameterization
9.	3D Coordinate Generation: Project optimized 2D molecular graphs into spatial, realistic 3D atomic coordinates.
10.	Bonding Patterns & Valence Verification: Audit structural hybridization (sp2 vs. sp3) following coordinate inflation to prevent geometric warping.
11.	Atomic Partial Charge Assignment: Compute robust quantum or semi-empirical electrostatic profiles via AM1-BCC or Gasteiger models over the explicit protonated systems.
12.	Forcefield Parameterization: Assign specific structural forcefield parameters (MMFF94, OPLS4, or GAFF2) to prepare for physical simulation.
13.	Energy Minimization: Apply gradient-descent algorithms to relieve structural strain and steer coordinates toward a true local potential energy minimum.
14.	Conformational Optimization / Ensemble Generation: Sample low-energy thermodynamic shapes for highly flexible marine leads and export the final ensemble as a multi-structure SDF or PDBQT file.


