Reproduction of [Lipid bilayers - Part II – Complex mixtures](https://cgmartini.nl/docs/tutorials/Martini3/LipidsII/) tutorial for Martini3

This tutorial covers:
1. Building and simulating a multi-component membrane without waiting for self-assembly
2. Adding a new lipid headgroup (lysyl-PG)

___
### Complex mixtures with Insane

`Insane` (INSert membraNE) can be found at this [repo](https://github.com/Tsjerk/Insane).

The custom `Insane` has settings for DBPC but not DLPC. We can supplement the database of lipid definitions that ships with `Insane` by creating a custom `lipid.dat` file with our definition for DLPC and referencing that using the `-dat` flag. Later, when we add a new lipid type, we will add that to our supplemental `lipids.dat` file. 

1. Create the file:
`touch lipids.dat`

```
[ diacylglycerol ]
;
; PROTOLIPID (diacylglycerol), 18 beads
;
; 1-3-4-6-7--9-10-11-12-13-14
;  \| |/  |
;   2 5   8-15-16-17-18-19-20
;
           0  .5   0   0  .5   0   0  .5   0   0   0   0   0   0   1   1   1   1   1   1
           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
          10   9   9   8   8   7   6   6   5   4   3   2   1   0   5   4   3   2   1   0
;          1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
M3.DLPC   0   -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   -
```
2. Build the starting membrane with `Insane`:
`insane -dat lipids.dat -salt 0.15 -x 20 -y 10 -z 10 -sol W -o membrane.gro -l DBPC:2 -l DLPC:1 -p topol.top`

This will create a `.top` file that needs the \#include statements to be modified into:
```
#include "../../martini_v300/martini_v3.0.0.itp" ; the particle definitions should be included first
#include "../../martini_v300/martini_v3.0.0_ffbonded_v2.itp" ; the general lipid bonded definitions should be included before the specific lipid topologies
#include "../../martini_v300/martini_v3.0.0_phospholipids_PC_v2.itp" ; remaining relevant .itps for molecules
#include "../../martini_v300/martini_v3.0.0_solvents_v1.itp" ; for water
#include "../../martini_v300/martini_v3.0.0_ions_v1.itp" ; for ions
```

3. Energy minimization:
```
gmx grompp -f em.mdp -c membrane.gro -p topol.top -o em.tpr
gmx mdrun -s em.tpr -v -c em.gro
```

4. Create an `index.ndx` file to group lipids together (as well as water an ions, solvent), to fit specified entries in the `md.mdp`
 file:

`gmx make_ndx -f em.gro`
```
2 | 3
4 | 5
name 6 Bilayer
name 7 Solvent
```

5. Quick production run:
```
gmx grompp -f md.mdp -p topol.top -c em.gro -n index.ndx -o phase_sep.tpr
gmx mdrun -v -deffnm phase_sep
```
___
### Analysis – Hexagonality

The *hexagonality* of lipid tail packing is a number that tells you how close the neighbors around each lipid are to forming a perfect hexagonal arrangement, here calculated along the phase transition, using only the first bead of each acyl tail. A tail counts as hexagonally packed if:

- its first bead has at least 6 neighbors within 6 Å, and
- all of its 6 closest neighbors also have at least 2 neighbors within 6 Å.

```
./hexagonality.py -f phase_sep.xtc -s phase_sep.tpr
qtgrace hexagonality.xvg
```
<img width="780" height="542" alt="Pasted image 20261028165441" src="https://github.com/user-attachments/assets/467ca678-d063-45db-9b57-9de41e4765c6" />

The y-axis is the fraction of tails that pass the aforementioned test (0 to 1). The system in this tutorial is 2:1 DBPC:DLPC, and DBPC is the gel-forming lipid. If every DBPC tail were hexagonally packed and no DLPC tail were, the fraction would be about 0.67. The system plateaus at ~0.61, expected of a well phase-separated system. We can estimate the temperature of each phase from the last frame of the simulation using the `kinetic-temperature.py` script, to check whether the gel (DBPC-rich) and liquid (DLPC-rich) regions both sit around 290 K.

`python3 kinetic-temperature.py'

```
DBPC    6300 beads  T = 287.5 K
DLPC    2688 beads  T = 296.0 K
W       8500 beads  T = 287.1 K
```
___
### Adding a new lipid headgroup: lysyl-PG

This section covers the parametrization of a new lipid headgroup, such as in cases where one isn't yet described in the Martini3 lipidome.

Lysyl-PG (lysyl-phosphatidylglycerol) is a membrane lipid in several gram-positive bacteria, such as Staphylococcus aureus, which confers a higher level of antimicrobial resistance to cationic host defense peptides[5]. Lysyl-PG consists of the ester of phosphatidylglycerol (PG) with the amino acid lysine.

<img width="555" height="241" alt="Pasted image 20261028181007" src="https://github.com/user-attachments/assets/4ddb544d-8aeb-4f72-b272-12bf2536af9d" />

#### Creating the lysyl-PG topology

```
mkdir new-headgroup
cd new-headgroup
touch lysilPG.top
```

Check the topology for POPG in `martini_v3.0.0_phospholipids_PG_v2.itp` and for Lys in `martini_v300/martini_v3.0.0_proteins/force_fields/martini3001/aminoacids.ff`. We'll also create a file named `lysilPG.top` and copy both topologies to it. The simplest way of doing this is by using the POPG topology as a template and adding the lysine parameters on top

```
[ moleculetype ]
; molname nrexcl
POLPG 1

[ atoms ]
; id type resnr residu atom cgnr charge
   1   P2    1  POLPG   BB   1    0   ;Lysine
   2  SC3    1  POLPG  SC1   2    0
   3 SQ4p    1  POLPG  SC2   3  1.0
   4  P4r    1  POLPG  GL0   4    0   ;POPG
   5   Q5    1  POLPG  PO4   5 -1.0
   6 SN4a    1  POLPG  GL1   6    0
   7 SN4a    1  POLPG  GL2   7    0
   8   C1    1  POLPG  C1A   8    0
(continues ...)
```

We'll also add it to the `lipids.dat` for Insane! and assemble the membrane:

```
[ lysylpg ]
;
; PROTOLIPID (lysyl-phosphatidylglycerol)
;
; 3  1--4--5--6--8--9--10--11--12--13     SC2  BB--GL0--PO4--GL1--C1A--(...)
;  \ |      \ |                              \ |           \  |
;   2        7--14--15--16--17--18--19       SC1             GL2--C1B--(...)
; 
               0  .5   0   0   0   0  .5   0   0   0   0   0   0   1   1   1   1   1   1  
               0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
               9   9  10   8   7   6   6   5   4   3   2   1   0   5   4   3   2   1   0
; bead number: 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
M3.POLPG   1   BB SC1 SC2 GL0 PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   -
; M3 indicates this template will be used when the Martini 3 FF is asked for (which is the default)
; The numeric '1' field is the lipid's overall charge
```

`insane -dat ../lipids.dat -salt 0.15 -l POPG:5 -l POLPG:1 -x 20 -y 10 -z 10 -sol W -charge auto -sol W -o membrane.gro -p topol.top`

Remember to modify the `topol.top` file to point at the right `.itp` files.

```
#include "../../../martini_v300/martini_v3.0.0.itp" ; the particle definitions should be included first
#include "../../../martini_v300/martini_v3.0.0_ffbonded_v2.itp" ; the general lipid bonded definitions should be included before the specific lipid topologies
#include "polpg.itp"
#include "../../../martini_v300/martini_v3.0.0_phospholipids_PC_v2.itp" ; remaining relevant .itps for molecules
#include "../../../martini_v300/martini_v3.0.0_solvents_v1.itp" ; for water
#include "../../../martini_v300/martini_v3.0.0_ions_v1.itp" ; for ions
```

Copy the `.mdp` files from before; run another energy minimization step and a production run (change temperature to 300 K and the total simulation time to 50 ns (nsteps = 2500000)).

`cd ../em.mdp ../md.mdp .`

Minimization:
```
gmx grompp -f em.mdp -c membrane.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

Make an `index` file again

Production run:
```
gmx grompp -f md.mdp -p topol.top -c em.gro -n index.ndx -o md.tpr
gmx 
```
