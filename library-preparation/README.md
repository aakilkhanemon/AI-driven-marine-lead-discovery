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

Table. Ligand preparation: all steps chronologically from input to docking-ready compound. Enlisted all Python libraries and modules.  

Phase / Step	Operation	Computational Model / Standard Reference	Python Library	Target Function / Modules	Library Expansion / Deduction
Phase I: 1	Deduplication	Canonical SMILES String Generation	rdkit.Chem	Chem.CanonSmiles, pandas.DataFrame.drop_duplicates	Deduction
(Removes redundant records)
Phase I: 2	Structure Sanitization	Standard Valence/Aromaticity Rules	rdkit.Chem	Chem.SanitizeMol, SanitizeFlags	Deduction
(Discards invalid/unfixable valencies)
Phase I: 3	Fragment Removal / Desalting	Ion-Stripping Filters (Na^+, Cl^-, etc.)
	rdkit.Chem.SaltRemover	SaltRemover.StripMol	Deduction
(Remove salts/solvents, reduces MW)
Phase I: 4	Functional Group Standardization	Functional Group Normalization Rules	rdkit.Chem.MolStandardize	rdkit.Chem.MolStandardize.rdMolStandardize	Neutral
(Standardizes tautomeric/functional forms)
Phase I: 5	Bioavailability Filtering	Lipinski’s Rule of 5 / Veber Criteria / REOS	rdkit.Chem.Descriptors	Descriptors.ExactMolWt, Descriptors.MolLogP, Lipinski.NumHDonors	Deduction
(Eliminates non-drug-like space)
Phase I: 6	PAINS Filtering	Pan-Assay Interference Compounds Filters	rdkit.Chem.FilterCatalog	FilterCatalogParams.FilterCategory.PAINS, FilterCatalog	Deduction
(Eliminates pan-assay interferers)
Phase II: 7	Physiological Protonation & Ionization	Physiological pH (7.4 \pm 0.2)
	Dimorphite-DL / Gypsum-DL	Dimorphite_DL.run_dimorphite_dl	Expansion
(Generates multiple protonation states at pH 7.4)
Phase II: 8	Stereochemical Fidelity & Enumeration	Chirality and Isomer Enumeration Engines	rdkit.Chem.EnumerateStereoisomers	EnumerateStereoisomers.EnumerateStereoisomers	Expansion
(Multiplies library via unassigned chiral centers)
Phase III: 9	3D Coordinate Generation	Distance Geometry / Knowledge-Based Embedding	rdkit.Chem.AllChem	AllChem.EmbedMolecule, AllChem.EmbedMultipleConfs	Neutral
(Transitions coordinates from 2D rightarrow 3D)
Phase III: 10	Bonding Patterns & Valence Verification	Valence / Hybridization Auditing	rdkit.Chem	Chem.DetectChemistryProblems, Mol.GetHybridization	Deduction
(Discards molecules warped during inflation)
Phase III: 11	Atomic Partial Charge Assignment	AM1-BCC / Gasteiger Quantum Models	openff.toolkit / rdkit.Chem	Molecule.assign_partial_charges(method='am1bcc')	Neutral
(Appends electronic parameters)
Phase III: 12	Forcefield Parameterization	MMFF94 / OPLS4 / GAFF2	openmmforcefields / AllChem	ForceFieldFactory, AllChem.MMFFGetMoleculeProperties	Neutral
(Appends physical parameter mechanics)
Phase III: 13	Energy Minimization	Gradient-Descent Algorithm / Conjugate Gradient	rdkit.Chem.AllChem	AllChem.MMFFOptimizeMolecule, AllChem.UFFOptimizeMolecule	Neutral
(Alters atomic coordinates to local minima)
Phase III: 14	Conformational Optimization	Low-Energy Ensemble Generation (SDF / PDBQT)	rdkit.Chem.AllChem / Omega	AllChem.PruneConformers, Chem.SDWriter	Expansion
(Generates low-energy multi-conformer ensembles)


