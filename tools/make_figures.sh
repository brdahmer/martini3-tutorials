#!/usr/bin/env bash
# Rebuild every figure made from an .xvg file. Run from anywhere:
#   bash tools/make_figures.sh
# The .xvg files stay where the tutorial commands write them; the PNGs go into
# each tutorial's figures/ folder.
set -euo pipefail
cd "$(dirname "$0")/.."
plot() { python3 tools/plot_xvg.py "$@"; }

RHO='Density (kg m$^{-3}$)'

# 00: DSPC self-assembly ------------------------------------------------------
T=00_dspc-spontaneous-assembly; F=$T/figures; mkdir -p "$F"
for d in x y z; do
  plot $T/dens_$d.xvg -o $F/density-$d.png --title "Density along $d after self-assembly" \
    --xlabel "Position along $d (nm)" --ylabel "$RHO"
  plot $T/dens_${d}_rot.xvg -o $F/density-$d-rotated.png --title "Density along $d after rotating the box" \
    --xlabel "Position along $d (nm)" --ylabel "$RHO"
done

E=$T/equilibration
# x and y are identical (semi-isotropic coupling, square box), so plot one line.
plot $E/box-xy.xvg -o $F/box-xy.png --ns --columns 1 --title "Box length in x and y during equilibration" \
  --subtitle "x and y are identical under semi-isotropic pressure coupling" --ylabel "Box length (nm)"
plot $E/area.xvg -o $F/area-per-lipid.png --ns --title "Area per lipid" \
  --xlabel "Time (ns)" --ylabel 'Area per lipid (nm$^2$)'
plot $E/errest.xvg -o $F/area-per-lipid-error.png --title "Error estimate of the area per lipid" \
  --subtitle 'Block averaging after 15 ns. Mean 0.5932 nm$^2$, error 0.0027 nm$^2$' \
  --xlabel "Block size (ps)" --ylabel 'Error estimate (nm$^2$)' --legend "Block averaging" "Analytical fit"
plot $E/p-density.xvg -o $F/phosphate-density.png --title "Phosphate (P*) density along z" \
  --xlabel "Position along z (nm)" --ylabel "$RHO"
plot $E/p-density-centered.xvg -o $F/phosphate-density-centered.png \
  --title "Phosphate (P*) density along z, bilayer centered" \
  --xlabel "Position along z (nm)" --ylabel "$RHO"
plot $E/msd.xvg -o $F/msd.png --ns --title "Lateral mean squared displacement of DSPC" \
  --subtitle 'D = 0.105 ± 0.019 × 10$^{-5}$ cm$^2$/s' --xlabel "Lag time τ (ns)" --ylabel 'MSD (nm$^2$)'

# 01: multicomponent membrane ------------------------------------------------
T=01_multicomponent-membrane; F=$T/figures; N=$T/new-headgroup; mkdir -p "$F"
HEX='Fraction of hexagonally packed tails'
plot $T/hexagonality.xvg -o $F/hexagonality.png --ns --title "Tail packing in DBPC:DLPC (2:1)" \
  --xlabel "Time (ns)" --ylabel "$HEX"
plot $N/hexagonality.xvg -o $F/new-headgroup-hexagonality.png --ns --title "Tail packing in POPG:POLPG (5:1)" \
  --xlabel "Time (ns)" --ylabel "$HEX"
plot $N/box-xy.xvg -o $F/new-headgroup-box-xy.png --ns --title "Box dimensions, POPG:POLPG (5:1)" \
  --ylabel "Box length (nm)" --legend "x" "y"
plot $N/p-density.xvg -o $F/new-headgroup-phosphate-density.png --title "Phosphate (P*) density along z, POPG:POLPG (5:1)" \
  --xlabel "Position along z (nm)" --ylabel "$RHO"
