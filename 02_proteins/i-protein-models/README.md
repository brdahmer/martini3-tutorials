# Protein structure models (CK1δ)

Reproduction of [Tutorial I.I: Setup of Martini 3 Protein Structure Models](https://cgmartini.nl/docs/tutorials/Martini3/ProteinsI/Tut1.html)

Covers the set up of different Martini 3 structure bias models for the protein kinase casein kinase delta 1 (CK1δ, PDB ID: 4JJR)

> Covers only the setup so far!

4JJR is crystallized as a dimer. While the official tutorial provides a clean initial structure, it also shows how to produce a valid `.pdb` file as follows:

1. Download the `.pdb` file:

`wget https://files.rcsb.org/download/4JJR.pdb`

2. Remove non-protein atoms:

`grep "^ATOM" 4JJR.pdb > 4JJR_clean.pdb`

3. Extract the monomer (removing chain B):

`grep " A " 4JJR_clean.pdb > 4JJR_clean_A.pdb`

4. The clean monomer file contains missing residues, which can be modelled with, for exemple, the [SWISS-MODEL webtool](https://swissmodel.expasy.org/) or [AlphaFold](https://alphafoldserver.com/). The output file from SWISS-MODEL is named `model_01.pdb`, renamed to `CK1d.pdb`.


5. Subsequently, check the protonation state with, for example, the [H++ webserver](http://newbiophysics.cs.vt.edu/H++/index.php) or [PROPKA](https://github.com/jensengroup/propka-3.0).

___
### Working Directory

For this tutorial, everything is done from one project folder. Each protein model we'll generate will be located in a folder within `models/`, because the output files of every model have the same names (`CK1d_cg.pdb`, `molecule_0.itp`, ...) and would overwrite each other. As the force field files are identical for all models, and for the sake of keeping things contained, the relevant files were moved to `ff/`.

```
i-protein-models/
├── atomistic/CK1d.pdb
├── ff/                           # see below for modifications to martini_v3.0.0.itp
│   ├── martini_v3.0.0.itp
│   ├── martini_v3.0.0_ions_v1.itp
│   └── martini_v3.0.0_solvents_v1.itp
├── mdp/
│   ├── minimization.mdp
│   ├── equilibration.mdp
│   └── dynamic.mdp
├── tools/OLIVES_v2.1.1b_M3.0.0.py
└── models/
    ├── UNBIASED/
    ├── EN/
    ├── GO/
    ├── OLIVES/
    └── OLIVES_DSSP/
```

Each model is generated from their own directory and the input is always `../../atomistic/CK1d.pdb`. For example:

```
mkdir -p models/<MODEL>
cd models/<MODEL>
```
___
### Force Field Files

1. Copy the three Martini3 files (here from `martini_v300`) into `ff/`, leaving the originals untouched:

```
mkdir ff
cp ../../martini_v300/martini_v3.0.0.itp ../../martini_v300/martini_v3.0.0_ions_v1.itp ../../martini_v300/martini_v3.0.0_solvents_v1.itp ff/
```

2. The GōMartini model needs the master `ff/martini_v3.0.0.itp` modified, but the modification is harmless for the remaining models. The additional atomtypes defined in the `go_atomtypes.itp` file which are produced during the GōMartini protein model generation (see later section) should be included at the end of the [ atomtypes ] directive as:

```
[ atomtypes ]
...
TX1er 36.0 0.000 A 0.0 0.0
W  72.0 0.000 A 0.0 0.0
SW 54.0 0.000 A 0.0 0.0
TW 36.0 0.000 A 0.0 0.0
U  24.0 0.000 A 0.0 0.0

#ifdef GO_VIRT
#include "go_atomtypes.itp"
#endif

[ nonbond_params ]
P6    P6  1 4.700000e-01    4.990000e+00
P6    P5  1 4.700000e-01    4.730000e+00
P6    P4  1 4.700000e-01    4.480000e+00
...
```

Similarly, the same needs to be done with `go_nbparams.itp` at the end of the [ nonbond_params ] directive (very end of the file):

```
...
TX2er  SQ1n  1 3.660000e-01    3.528000e+00
TX2er  TQ1n  1 3.520000e-01    5.158000e+00
TX1er   Q1n  1 3.950000e-01    1.981000e+00
TX1er  SQ1n  1 3.780000e-01    3.098000e+00
TX1er  TQ1n  1 3.660000e-01    4.422000e+00

#ifdef GO_VIRT
#include "go_nbparams.itp"
#endif
```
___
### Unbiased Martini 3 Protein Model

See https://vermouth-martinize.readthedocs.io/en/v0.14.0/tutorials/basic_usage.html

This model is for illustration only and is the starting point of OLIVES + DSSP (see later section).

1. Using the [martinize2](https://github.com/marrink-lab/vermouth-martinize) program to generate GROMACS-compatible Coarse Grained structure and topology files:

```
mkdir -p models/UNBIASED
cd models/UNBIASED
martinize2 -f ../../atomistic/CK1d.pdb -x CK1d_cg.pdb -o CK1d_only.top -ff martini3001 -p backbone -dssp
```

where: 
- `-ff martini3001` specifies the forcefield
- `-p backbone` generates position restraints on the protein backbone
- `-dssp` adds secondary structure-dependent bonded potentials (if left blank, uses `mdtraj` for this)
this generates, in `models/UNBIASED/`:
- a CG protein structure (`CK1d_cg.pdb`) 
- a GROMACS topology file (`CK1d_only.top`) 
- a protein .itp file (`molecule_0.itp`)

2. Continue with *Solvation, Topology and Simulation* below.

___

### Elastic Network Model

See https://vermouth-martinize.readthedocs.io/en/v0.14.0/tutorials/elastic_networks.html

1. Is generated again with `martinize2`, using the `-elastic` flag with additional parameters (`-ef`, for force constant of harmonic potentials \[kept above $700 kJ/(mol\ nm^2)$] `-el`, lower cutoff for generation of elastic bonds, `-eu`, upper cutoff)

```
mkdir -p models/EN
cd models/EN
martinize2 -f ../../atomistic/CK1d.pdb -x CK1d_cg.pdb -o CK1d_only.top -ff martini3001 -p backbone -dssp -elastic -el 0 -eu 0.85
```

Same output files as the unbiased model (`CK1d_cg.pdb`, `CK1d_only.top`, `molecule_0.itp`), now with the elastic bonds in the `.itp`.

2. Continue with *Solvation, Topology and Simulation* below.
___

### GōMartini Model

See https://vermouth-martinize.readthedocs.io/en/v0.14.0/tutorials/go_models.html

The generation of a Gō-like network requires a contact map based on the atomistic input structure. Martinize2 can calculate this automatically by giving the `-go` flag without any argument.

1. Then, again with martinize2:

```
mkdir -p models/GO
cd models/GO
martinize2 -f ../../atomistic/CK1d.pdb -x CK1d_cg.pdb -o CK1d_only.top -ff martini3001 -p backbone -dssp -go -go-low 0.3 -go-up 1.1 -go-eps 9.414
```

Note: without `-go` the other `-go-*` flags don't switch the model on. With `-go` the protein itp is named `molecule.itp` (not `molecule_0.itp`), which matters for the topology later. If the contact map comes from the Gō web server instead, give it as `-go contact_map.out`.

Along with the standard files as before, the generation of a Gō-like network via martinize2 also produces the files `go_atomtypes.itp` and `go_nbparams.itp`, which define  the new virtual sites as atoms for your system and the specific non-bonded interactions between them, respectively. They stay in `models/GO/` (keep filenames as is).

2. In the `.top` file of your working system, the include statement `#define GO_VIRT` needs to be present (at the top, before the `#include` of `martini_v3.0.0.itp`). It is written in the topology step below.

3. Continue with *Solvation, Topology and Simulation* below.
___
### OLIVES

The OLIVES protein model is usually employed without the `-dssp` restraints to enable conformational changes and partial unfolding of domains (using the `OLIVES_v2.1.1b_M3.0.0.py` from the [OLIVES](https://github.com/Martini-Force-Field-Initiative/OLIVES) GitHub, kept in `tools/`). Two versions are generated:

1. OLIVES (no DSSP), which needs its own `martinize2` run without `-dssp`. Recent versions of `martinize2` refuse to run without `-dssp` or `-ss`, so give an explicit coil with `-ss C`:

```
mkdir -p models/OLIVES
cd models/OLIVES
martinize2 -f ../../atomistic/CK1d.pdb -x CK1d_cg.pdb -o CK1d_only.top -ff martini3001 -p backbone -ss C
python3 ../../tools/OLIVES_v2.1.1b_M3.0.0.py -c "CK1d_cg.pdb" -i "molecule_0.itp"
```

2. OLIVES + DSSP, same as above but with the `-dssp` flag (the `martinize2` command is the same as the unbiased model):

```
mkdir -p models/OLIVES_DSSP
cd models/OLIVES_DSSP
martinize2 -f ../../atomistic/CK1d.pdb -x CK1d_cg.pdb -o CK1d_only.top -ff martini3001 -p backbone -dssp
python3 ../../tools/OLIVES_v2.1.1b_M3.0.0.py -c "CK1d_cg.pdb" -i "molecule_0.itp"
```

The OLIVES script modifies the `.itp` it is given, and needs `mdtraj` and `pandas`. Comparing `models/OLIVES_DSSP/molecule_0.itp` with `models/UNBIASED/molecule_0.itp` shows what OLIVES added.

3. Continue with *Solvation, Topology and Simulation* below, in each OLIVES folder.

____

### Solvation, Topology and Simulation

Done for each model, from inside its own folder (`cd models/<MODEL>`). What changes between models:

| model | protein itp | `#define GO_VIRT` | `insane -charge` | extra `mdrun` flags |
|---|---|---|---|---|
| UNBIASED | `molecule_0.itp` | no | not needed | none |
| EN | `molecule_0.itp` | no | not needed | none |
| GO | `molecule.itp` | yes | `14` | none |
| OLIVES | `molecule_0.itp` | no | not needed | `-noddcheck -rdd 2.0` |
| OLIVES_DSSP | `molecule_0.itp` | no | not needed | `-noddcheck -rdd 2.0` |

1. Net charge of the CG protein (sum of the charges in `[ atoms ]`, 14 here, which also matches H++):

```
awk '/^\[/{sec=$2} sec=="atoms" && NF>=7 && $1!~/^;/{q+=$7} END{print q}' molecule_0.itp
```

2. Box, water and 0.15 M NaCl with `insane`:

```
insane -f CK1d_cg.pdb -o CG.gro -p insane.top -pbc cubic -box 10,10,10 -salt 0.15 -sol W -d 0
```

   - For GōMartini add `-charge 14` (the real net charge). `insane` counts the virtual sites as residues and reports double the charge (28), which gives the wrong number of ions.
   - For the other models `insane` finds the charge correctly (14).

3. Martini 3 ions have no charge sign, so remove it from the `.gro` and the topology (width-preserving, macOS `sed`):

```
sed -i '' 's/NA+/NA /g; s/CL-/CL /g' CG.gro
sed -i '' 's/^NA+/NA/; s/^CL-/CL/' insane.top
```

4. Master topology, `system.top`, written by hand: the includes at the top, then the `[ system ]` and `[ molecules ]` sections with the counts from `insane.top` (example for EN):

```
#include "../../ff/martini_v3.0.0.itp"
#include "../../ff/martini_v3.0.0_ions_v1.itp"
#include "../../ff/martini_v3.0.0_solvents_v1.itp"
#include "./molecule_0.itp"

[ system ]
; name
CK1d

[ molecules ]
; name  number
molecule_0 1
W   8203
NA    83
CL    97
```

   - For GōMartini add `#define GO_VIRT` as the first line, and use `./molecule.itp` and `molecule 1`.
   - The molecule name in `[ molecules ]` has to match the `[ moleculetype ]` name inside the protein `.itp` (`molecule_0`, or `molecule` for Gō), not the `Protein` that `insane` writes.
   - The ion balance (CL - NA) should equal the protein charge, 14, in every model. In GōMartini the counts are 8207 W, 83 NA and 97 CL.

5. Minimization, then equilibration, then production. The `.mdp` files are in `../../mdp/`; `-I .` makes `grompp` look in the model folder (where the Gō `.itp` files are), so run everything from there. Add the extra `mdrun` flags from the table:

```
gmx grompp -p system.top -c CG.gro -r CG.gro -f ../../mdp/minimization.mdp -o minimization.tpr
gmx mdrun -deffnm minimization -v

gmx grompp -p system.top -c minimization.gro -r CG.gro -f ../../mdp/equilibration.mdp -o equilibration.tpr
gmx mdrun -deffnm equilibration -v

gmx grompp -p system.top -c equilibration.gro -f ../../mdp/dynamic.mdp -o dynamic.tpr -maxwarn 1
gmx mdrun -deffnm dynamic -v
```

For replicas (the tutorial recommends at least 3 of 250 ns), repeat the last two commands with other names (`dynamic_2`, `dynamic_3`); `dynamic.mdp` regenerates velocities with a random seed.
___
### MDP Files

The tutorial's link only provides the generic Martini 3 `.mdp` (a membrane-oriented template), so the three files in `mdp/` are derived from it:

- `minimization.mdp`: steepest descent, up to 10000 steps
- `equilibration.mdp`: 1 ns, `define = -DPOSRES`, `dt = 0.02`, `gen_vel = yes`
- `dynamic.mdp`: 250 ns (12500000 steps), no restraints, `gen_vel = yes` with a random seed

and the changes from the template are:
- `Pcoupltype = isotropic` with a single `ref_p` and compressibility (the template's `semiisotropic` is for bilayers), `c-rescale`, `tau_p = 4.0`
- `tc-grps = Protein non-Protein` (default index groups, instead of `Solvent Rest`, which needs an `.ndx` file), `comm-grps` removed
- `ref_t` / `gen_temp` = 310 K as in the template, adjust if needed
- the obsolete `ns_type` and `energygrps` removed
