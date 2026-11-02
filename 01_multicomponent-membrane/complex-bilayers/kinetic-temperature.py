import MDAnalysis as mda

u = mda.Universe("phase_sep.tpr", "phase_sep.gro")
kB = 0.0083144626  # kJ/mol/K

def temp(ag):
    v = ag.velocities / 10.0            # Å/ps -> nm/ps
    m = ag.masses[:, None]              # amu
    vcom = (m * v).sum(0) / m.sum()
    ke = 0.5 * (m * (v - vcom) ** 2).sum()   # amu*nm²/ps² = kJ/mol
    return 2 * ke / ((3 * ag.n_atoms - 3) * kB)

for name in ("DBPC", "DLPC", "W"):
    ag = u.select_atoms(f"resname {name}")
    print(f"{name:5s} {ag.n_atoms:6d} beads  T = {temp(ag):.1f} K")