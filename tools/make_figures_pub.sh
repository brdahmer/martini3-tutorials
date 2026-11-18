#!/usr/bin/env bash
# Publication-style versions of the figures made by make_figures.sh.
# Output goes to each tutorial's figures-publication/ folder (PNG + PDF), so the
# originals in figures/ stay as they are for comparison. Run from anywhere:
#   bash tools/make_figures_pub.sh
set -euo pipefail
cd "$(dirname "$0")/.."
plot() { python3 tools/plot_xvg_pub.py "$@"; }

RHO='Density (kg m$^{-3}$)'

# 00: DSPC self-assembly ------------------------------------------------------
T=00_dspc-spontaneous-assembly; F=$T/figures-publication
for d in x y z; do
  plot $T/dens_$d.xvg -o $F/density-$d --xlabel "$d (nm)" --ylabel "$RHO"
  plot $T/dens_${d}_rot.xvg -o $F/density-$d-rotated --xlabel "$d (nm)" --ylabel "$RHO"
done

E=$T/equilibration
plot $E/box-xy.xvg -o $F/box-xy --ns --columns 1 --xlabel "Time (ns)" --ylabel "Box length x = y (nm)"
plot $E/area.xvg -o $F/area-per-lipid --ns --xlabel "Time (ns)" --ylabel 'Area per lipid (nm$^2$)'
plot $E/errest.xvg -o $F/area-per-lipid-error --ns --xlabel "Block size (ns)" --ylabel 'Error estimate (nm$^2$)' \
  --legend "Block averaging" "Analytical fit" --note 'APL = 0.5932 ± 0.0027 nm$^2$'
plot $E/p-density.xvg -o $F/phosphate-density --xlabel "z (nm)" --ylabel "$RHO"
plot $E/p-density-centered.xvg -o $F/phosphate-density-centered --xlabel "z, bilayer centered (nm)" --ylabel "$RHO"
plot $E/msd.xvg -o $F/msd --ns --xlabel "Lag time (ns)" --ylabel 'Lateral MSD (nm$^2$)' \
  --note 'D = (0.105 ± 0.019) × 10$^{-5}$ cm$^2$ s$^{-1}$'

# 01: multicomponent membrane ------------------------------------------------
T=01_multicomponent-membrane; F=$T/figures-publication; N=$T/new-headgroup
HEX='Hexagonally packed tails (fraction)'
plot $T/hexagonality.xvg -o $F/hexagonality --ns --xlabel "Time (ns)" --ylabel "$HEX" --note "DBPC:DLPC 2:1"
plot $N/hexagonality.xvg -o $F/new-headgroup-hexagonality --ns --xlabel "Time (ns)" --ylabel "$HEX" --note "POPG:POLPG 5:1"
plot $N/box-xy.xvg -o $F/new-headgroup-box-xy --ns --xlabel "Time (ns)" --ylabel "Box length (nm)" --legend "x" "y"
plot $N/p-density.xvg -o $F/new-headgroup-phosphate-density --xlabel "z (nm)" --ylabel "$RHO"

# 02: transmembrane peptide tilt ---------------------------------------------
(cd 02_proteins/ii-transmembrane-proteins/analysis && python3 plot_tilt_pub.py)
