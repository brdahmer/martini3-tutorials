import COBY
COBY.COBY(
### Def. box size (xyz, nm) and type.
box = [7, 7, 9],
box_type = "rectangular",
### 'membrane' creates a simple POPC membrane
membrane = "lipid:POPC",
### 'protein' inserts a protein at the center of the system
protein = "file:kalp.pdb moleculetypes:molecule_0",
### 'solvation' solvates the system using the default settings [water] 55.56 M & [NaCl] = .15 M
solvation = "default",
### Include .itps for charge calcs and to write .top file.
itp_input = [
"include:../../../martini_v300/martini_v3.0.0.itp",
"include:../../../martini_v300/martini_v3.0.0_ions_v1.itp",
"include:../../../martini_v300/martini_v3.0.0_phospholipids_v1.itp",
"include:../../../martini_v300/martini_v3.0.0_solvents_v1.itp",
"include:../peptide-martinize/molecule_0.itp", ],
### Designates the system name
sn = "SimpleMembrane",
### File writing
out_sys = "system.gro", out_top = "system.top", out_log = "log.log")
