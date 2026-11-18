#!/usr/bin/env bash
# setup_ck1d_models.sh
#
# Builds one self-contained folder per Martini 3 structure bias model
# (unbiased, EN, GoMartini, OLIVES, OLIVES+DSSP) from a single atomistic CK1d.pdb,
# following Martini Tutorial I.I. Written to work with macOS (bash 3.2, BSD sed).
#
# Expected layout before running:
#   $ROOT/atomistic/CK1d.pdb            (model_01.pdb, renamed)
#   $FF_DIR/martini_v3.0.0{,_ions_v1,_solvents_v1}.itp   (source, never modified; default $ROOT/ff_source)
#   -> copied once to $FF_SHARED (default $ROOT/ff) and patched with the GO_VIRT includes
#   $ROOT/tools/OLIVES_v2.1.1b_M3.0.0.py   (only for the OLIVES models)
#   $ROOT/contact_map.out                  (only for GoMartini; see below)
#   $ROOT/mdp/{minimization,equilibration,dynamic}.mdp   (only if RUN_MD=1)
#
# Usage:
#   ./setup_ck1d_models.sh                 # build all five models
#   ./setup_ck1d_models.sh EN GO           # build a subset
#   RUN_MD=1 ./setup_ck1d_models.sh EN     # also minimize/equilibrate/produce
#
# Model names: UNBIASED  EN  GO  OLIVES  OLIVES_DSSP  (all built by default)

set -euo pipefail

ROOT="${ROOT:-$PWD}"
PDB="${PDB:-$ROOT/atomistic/CK1d.pdb}"
FF_DIR="${FF_DIR:-$ROOT/ff_source}"     # original martini_v3.0.0*.itp files (never modified)
FF_SHARED="${FF_SHARED:-$ROOT/ff}"      # one patched copy, shared by all models
MDP_DIR="${MDP_DIR:-$ROOT/mdp}"
OUT_DIR="${OUT_DIR:-$ROOT/models}"
CONTACT_MAP="${CONTACT_MAP:-$ROOT/contact_map.out}"
OLIVES_SCRIPT="${OLIVES_SCRIPT:-$ROOT/tools/OLIVES_v2.1.1b_M3.0.0.py}"
GMX="${GMX:-gmx}"
RUN_MD="${RUN_MD:-0}"

# Model parameters from the tutorial
EN_LOWER="${EN_LOWER:-0}"
EN_UPPER="${EN_UPPER:-0.85}"
GO_LOW="${GO_LOW:-0.3}"
GO_UP="${GO_UP:-1.1}"
GO_EPS="${GO_EPS:-9.414}"

# insane charge override for GO only. Leave empty to use the protein charge computed
# from the itp (insane double counts the Go virtual sites and would report 2x).
INSANE_CHARGE_GO="${INSANE_CHARGE_GO:-}"

BOX="${BOX:-10,10,10}"
SALT="${SALT:-0.15}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo ">>> $*"; }

need() { command -v "$1" >/dev/null 2>&1 || die "'$1' not found in PATH"; }

need martinize2
need insane
need python3
[ "$RUN_MD" = "1" ] && need "$GMX"

[ -f "$PDB" ] || die "atomistic input not found: $PDB"
# A force field file must exist either in the shared folder (already set up) or in the source folder
for f in martini_v3.0.0.itp martini_v3.0.0_ions_v1.itp martini_v3.0.0_solvents_v1.itp; do
  [ -f "$FF_SHARED/$f" ] || [ -f "$FF_DIR/$f" ] \
    || die "missing force field file $f (looked in $FF_SHARED and $FF_DIR; set FF_DIR to the folder with the originals)"
done

if [ "$#" -gt 0 ]; then
  MODELS="$*"
else
  MODELS="UNBIASED EN GO OLIVES OLIVES_DSSP"
fi

# ---------------------------------------------------------------------------
# Step 1: martinize2 (+ OLIVES / Go handling) inside the model folder
# ---------------------------------------------------------------------------
coarse_grain() (
  model=$1
  dir=$2
  cd "$dir"

  base="-f $PDB -x CK1d_cg.pdb -o CK1d_only.top -ff martini3001 -p backbone"

  case "$model" in
    UNBIASED)     flags="$base -dssp" ;;
    EN)           flags="$base -dssp -elastic -el $EN_LOWER -eu $EN_UPPER" ;;
    OLIVES)       flags="$base -ss C" ;;   # no DSSP; newer martinize2 requires an explicit choice
    OLIVES_DSSP)  flags="$base -dssp" ;;
    GO)
      if [ -f "$CONTACT_MAP" ]; then
        gomap="$CONTACT_MAP"
      else
        gomap=""
        log "No contact map at $CONTACT_MAP; trying bare -go (needs martinize2 >= 0.13.0)"
      fi
      flags="$base -dssp -go $gomap -go-low $GO_LOW -go-up $GO_UP -go-eps $GO_EPS"
      ;;
    *) die "unknown model: $model" ;;
  esac

  # Remove outputs of earlier runs so the itp detection below is unambiguous
  rm -f molecule*.itp go_atomtypes.itp go_nbparams.itp

  log "[$model] martinize2"
  # shellcheck disable=SC2086
  martinize2 $flags 2>&1 | tee martinize.log

  # martinize2 names the protein itp molecule_0.itp normally, but molecule.itp with -go
  itp=""
  for c in molecule_0.itp molecule.itp; do
    if [ -f "$c" ]; then itp="$c"; break; fi
  done
  [ -n "$itp" ] || die "[$model] martinize2 produced no molecule*.itp (see martinize.log)"
  echo "$itp" > .itp_name
  log "[$model] protein itp: $itp"

  case "$model" in
    OLIVES|OLIVES_DSSP)
      [ -f "$OLIVES_SCRIPT" ] || die "OLIVES script not found: $OLIVES_SCRIPT"
      log "[$model] OLIVES"
      python3 "$OLIVES_SCRIPT" -c CK1d_cg.pdb -i "$itp" 2>&1 | tee olives.log
      ;;
  esac
)

# ---------------------------------------------------------------------------
# Step 2: copy force field, apply GoMartini edits to the LOCAL copy only
# ---------------------------------------------------------------------------
prepare_shared_ff() (
  mkdir -p "$FF_SHARED"
  for f in martini_v3.0.0.itp martini_v3.0.0_ions_v1.itp martini_v3.0.0_solvents_v1.itp; do
    # never overwrite an existing (possibly already patched) shared copy
    [ -f "$FF_SHARED/$f" ] || cp "$FF_DIR/$f" "$FF_SHARED/$f"
  done

  # Add the GO_VIRT includes to the shared master itp. The includes are guarded by
  # #ifdef GO_VIRT, so they do nothing for models whose system.top does not define it.
  # The bare file names are resolved from the directory grompp is run in (the model
  # folder), where martinize2 wrote go_atomtypes.itp / go_nbparams.itp.
  # Python instead of sed: portable between macOS and Linux, and idempotent.
  python3 - "$FF_SHARED/martini_v3.0.0.itp" <<'EOF'
import sys
p = sys.argv[1]
s = open(p).read()
if 'go_atomtypes.itp' in s:
    sys.exit(0)
marker = '[ nonbond_params ]'
if marker not in s:
    sys.exit('could not find [ nonbond_params ] in ' + p)
s = s.replace(marker,
              '#ifdef GO_VIRT\n#include "go_atomtypes.itp"\n#endif\n\n' + marker, 1)
s += '\n#ifdef GO_VIRT\n#include "go_nbparams.itp"\n#endif\n'
open(p, 'w').write(s)
EOF
)

prepare_ff() (
  model=$1
  dir=$2
  if [ "$model" = "GO" ]; then
    [ -f "$dir/go_atomtypes.itp" ] && [ -f "$dir/go_nbparams.itp" ] \
      || die "[GO] go_atomtypes.itp / go_nbparams.itp missing after martinize2"
  fi
)

# ---------------------------------------------------------------------------
# Step 3: insane box + master topology
# ---------------------------------------------------------------------------
build_system() (
  model=$1
  dir=$2
  cd "$dir"
  itp=$(cat .itp_name)

  # Net charge of the CG protein, summed from the [ atoms ] section of the itp
  netq=$(awk '
    /^\[/ { sec = $2 }
    sec == "atoms" && NF >= 7 && $1 !~ /^;/ { q += $7 }
    END { print q + 0 }' "$itp")

  # Go virtual sites make insane count each residue twice (it reports 2x the charge),
  # so for GO we pass the true charge explicitly.
  charge_flag=""
  if [ "$model" = "GO" ]; then
    charge_flag="-charge ${INSANE_CHARGE_GO:-$netq}"
  fi

  log "[$model] insane"
  # shellcheck disable=SC2086
  insane -f CK1d_cg.pdb -o CG.gro -p insane.top -pbc cubic -box "$BOX" \
         -salt "$SALT" -sol W -d 0 $charge_flag 2>&1 | tee insane.log

  # Martini 3 ion names have no charge sign (NA, CL, not NA+, CL-).
  # Same-width replacement keeps the fixed-column .gro valid.
  sed -e 's/NA+/NA /g' -e 's/CL-/CL /g' CG.gro > CG.tmp && mv CG.tmp CG.gro

  # Molecule name must match the [ moleculetype ] name inside the protein itp.
  mol=$(awk '/^\[ *moleculetype *\]/{f=1;next} f && NF && $1 !~ /^;/ {print $1; exit}' "$itp")
  [ -n "$mol" ] || die "[$model] could not read moleculetype name from $itp"

  # Path from the model folder to the shared force field folder
  ffrel=$(python3 -c 'import os,sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))' "$FF_SHARED" "$PWD")

  {
    [ "$model" = "GO" ] && echo "#define GO_VIRT"
    echo "#include \"$ffrel/martini_v3.0.0.itp\""
    echo "#include \"$ffrel/martini_v3.0.0_ions_v1.itp\""
    echo "#include \"$ffrel/martini_v3.0.0_solvents_v1.itp\""
    echo "#include \"./$itp\""
    echo
    awk -v mol="$mol" '
      /^#include/ { next }
      /^\[ *molecules *\]/ { inmol = 1; print; next }
      inmol && NF && $1 !~ /^;/ && !done { $1 = mol; done = 1 }
      { print }' insane.top \
      | sed -e 's/^NA+/NA/' -e 's/^CL-/CL/'
  } > system.top

  # Sanity check: net charge of the CG protein, and ion counts
  log "[$model] protein net charge from $itp: $netq"
  log "[$model] [ molecules ] section:"
  awk '/^\[ *molecules *\]/{f=1} f' system.top
)

# ---------------------------------------------------------------------------
# Step 4 (optional): minimization -> equilibration -> production
# ---------------------------------------------------------------------------
run_md() (
  model=$1
  dir=$2
  cd "$dir"
  for f in minimization equilibration dynamic; do
    [ -f "$MDP_DIR/$f.mdp" ] || die "missing $MDP_DIR/$f.mdp"
  done

  extra=""
  case "$model" in
    OLIVES|OLIVES_DSSP) extra="-noddcheck -rdd 2.0" ;;
  esac

  log "[$model] minimization"
  $GMX grompp -I . -p system.top -c CG.gro -r CG.gro -f "$MDP_DIR/minimization.mdp" -o minimization.tpr
  # shellcheck disable=SC2086
  $GMX mdrun -deffnm minimization -v $extra

  log "[$model] equilibration"
  $GMX grompp -I . -p system.top -c minimization.gro -r CG.gro -f "$MDP_DIR/equilibration.mdp" -o equilibration.tpr
  # shellcheck disable=SC2086
  $GMX mdrun -deffnm equilibration -v $extra

  log "[$model] production"
  $GMX grompp -I . -p system.top -c equilibration.gro -f "$MDP_DIR/dynamic.mdp" -o dynamic.tpr -maxwarn 1
  # shellcheck disable=SC2086
  $GMX mdrun -deffnm dynamic -v $extra
)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
mkdir -p "$OUT_DIR"
prepare_shared_ff

for model in $MODELS; do
  dir="$OUT_DIR/$model"
  mkdir -p "$dir"
  log "===== $model -> $dir ====="
  coarse_grain "$model" "$dir"
  prepare_ff   "$model" "$dir"
  build_system "$model" "$dir"
  if [ "$RUN_MD" = "1" ]; then
    run_md "$model" "$dir"
  fi
done

log "Done. Models built: $MODELS"