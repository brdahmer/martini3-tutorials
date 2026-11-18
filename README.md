# Martini 3 tutorials

Worked reproductions of the official [Martini 3 tutorials](https://cgmartini.nl/docs/tutorials/Martini3/) on lipids and proteins for coarse-grained molecular dynamics with GROMACS.

## Contents

| Folder | Tutorial | What it covers |
|---|---|---|
| [`00_dspc-spontaneous-assembly/`](00_dspc-spontaneous-assembly/) | [Lipid bilayers I – Self-assembly](https://cgmartini.nl/docs/tutorials/Martini3/LipidsI/) | Spontaneous self-assembly of a 128-lipid DSPC bilayer; equilibration; density profiles, area per lipid, lateral diffusion (MSD) and tail order parameters |
| [`01_multicomponent-membrane/`](01_multicomponent-membrane/) | [Lipid bilayers II – Complex mixtures](https://cgmartini.nl/docs/tutorials/Martini3/LipidsII/) | Building a DBPC/DLPC mixed membrane with `insane`, phase separation, hexagonal packing analysis, and parametrising a new lipid headgroup (lysyl-PG) |
| [`02_proteins/i-protein-models/`](02_proteins/i-protein-models/) | [Proteins I – Structure models](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html) | Coarse-graining casein kinase 1δ (PDB 4JJR) with `martinize2` and setting up different Martini 3 structure-bias models |
| [`02_proteins/ii-transmembrane-proteins/`](02_proteins/ii-transmembrane-proteins/) | [Proteins IIa – Transmembrane peptides](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsIIa/) | A KALP peptide in a POPC bilayer: water biasing, building the membrane four ways (self-assembly, `insane`, COBY, TS2CG), helix tilt analysis, a thinner DLPC bilayer and a dimerization run |
| [`martini_v300/`](martini_v300/) | – | Martini 3.0.0 force-field parameters used by all tutorials (third-party, see below) |


## Force field and third-party files

The files in `martini_v300/` are the [Martini 3 release](https://cgmartini.nl) parameters together with the [Martini 3 lipid parameters](https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters) in one folder for convenience.

- Souza, P. C. T. *et al.* Martini 3: a general purpose force field for coarse-grained molecular dynamics. *Nat. Methods* **18**, 382–388 (2021). [doi:10.1038/s41592-021-01098-3](https://doi.org/10.1038/s41592-021-01098-3)
- Wassenaar, T. A. *et al.* Computational lipidomics with *insane*: a versatile tool for generating custom membranes for molecular simulations. *J. Chem. Theory Comput.* **11**, 2144–2155 (2015). [doi:10.1021/acs.jctc.5b00209](https://doi.org/10.1021/acs.jctc.5b00209)
