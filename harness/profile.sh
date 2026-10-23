#!/usr/bin/env bash
# Usage: profile.sh baseline|optimized [--only <cell>] [--dry-run]
#   Profile VarDictRust (Rust engine only, no oracle) on a fixed set of named cells with the three profiler builds
#   of crates/vardictrust (features profile-cpu / profile-heap / profile-dhat) and rank the result:
#     scratch/profiles/<tag>/<cell>/{cpu-*,heap*,dhat-*,time.txt,cmd.txt,out-*}  +  scratch/profiles/<tag>/catalog.md
#   (catalog = harness/profile_catalog.py; SVG flamegraphs per cell are linked from it.)
#
# Cells (commands copied from existing harness files, never invented):
#   tiny s1a            germline slices, harness/slices.tsv rows (run_slice.sh columns)   [dhat]
#   s8 s9               indel/SV slices, harness/slices.tsv                               [dhat]
#   s3                  somatic tumor|normal pair, harness/slices.tsv                     [dhat]
#   rna-i               extended m1: -R chr8:90534611-90620070 -f 0.1 -i (run_extended_coverage.sh)  [dhat]
#   hg002               HG002 chr20:10M-30M slice (run_hg002.sh arguments)                [no dhat]
#   hg008               HG008 T|N slice, engine.bed (benchmarks/hg008-v04-guide.md cmd)   [no dhat]
#   wgs-scale           scratch/g4/seqc2/HCC1395.bam, bounded 10 Mb -R window at -th 8 (scale_bench.sh shape)  [no dhat]
#   cells whose inputs are missing are SKIPped (an error if named via --only). dhat only ever runs on slice cells
#   (its live-block map costs GBs: the systemd-oomd failure mode); HG002/HG008/WGS get cpu + heap only.
#
# Safety: every real run goes through harness/courtesy.sh (nice 19, ionice idle, disk/RAM/maintenance-window gate;
# a refusal aborts the whole run, never bypassed) and, for dhat and for hg002/hg008/wgs-scale, inside
# `systemd-run --user --scope -p MemoryMax=14G -p MemorySwapMax=0` when available. Threads are capped at 8
# (Host Policy A3): PROFILE_THREADS (default 8). Each run appends `/usr/bin/time -v` to <cell>/time.txt.
# Refuses to write under benchmarks/evidence/.
#
# Binaries: PROFILE_BIN_CPU / PROFILE_BIN_HEAP / PROFILE_BIN_DHAT (a single PROFILE_BIN overrides all three; used by
# the stub tests). Default: <base>/profile-<cpu|heap|dhat>/profiling/vardictrust with base = PROFILE_TARGET_BASE, else
# $CARGO_TARGET_DIR, else <repo>/target (so the default build trees sit under the git-ignored /target and a
# `git add -A` can never stage them; .gitignore also carries /target-*/); `cargo build --profile profiling -p vardictrust --features profile-<x>`
# is run for each default binary (incremental: a no-op when current, so a stale binary is never profiled).
# --dry-run prints the plan only (no build, no courtesy gate, nothing written).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
ROOT_DIR=$PWD

usage() { echo "usage: profile.sh baseline|optimized [--only <cell>] [--dry-run]" >&2; exit 2; }
TAG=${1:-}; [ $# -gt 0 ] && shift || true
case "$TAG" in baseline|optimized) ;; *) usage ;; esac
ONLY=""; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --only) [ $# -ge 2 ] || usage; ONLY=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    *) usage ;;
  esac
done

ALL_CELLS="tiny s1a s8 s9 s3 rna-i hg002 hg008 wgs-scale"
if [ -n "$ONLY" ]; then
  case " $ALL_CELLS " in *" $ONLY "*) ;; *) echo "profile: unknown cell '$ONLY' (known: $ALL_CELLS)" >&2; exit 2 ;; esac
fi

TH=${PROFILE_THREADS:-8}
[[ "$TH" =~ ^[0-9]+$ ]] && [ "$TH" -ge 1 ] && [ "$TH" -le 8 ] \
  || { echo "profile: PROFILE_THREADS=$TH out of range 1..8 (Host Policy A3 cap)" >&2; exit 2; }

OUT_ROOT=${PROFILE_OUT:-scratch/profiles}
OUT="$OUT_ROOT/$TAG"
EVIDENCE=$(realpath -m benchmarks/evidence)
case "$(realpath -m "$OUT")/" in
  "$EVIDENCE"/*) echo "profile: refusing to write under benchmarks/evidence/ ($OUT)" >&2; exit 2 ;;
esac

# --- cells -----------------------------------------------------------------------------------------------
# cell_spec <name> sets C_REF, C_BAM (a|b = tumor|normal), C_PRE (args before -th), C_POST (BED positional),
# C_DHAT (0/1), C_BIG (0/1), C_EST (courtesy EST_MIN minutes), C_INPUTS (files that must exist).
slice_spec() {   # same columns run_slice.sh reads from harness/slices.tsv; compare.sh's Rust argv shape
  local row EXTRA SAMPLE BED; row=$(awk -F'\t' -v id="$1" '$1==id' harness/slices.tsv)
  [ -n "$row" ] || { echo "profile: no harness/slices.tsv row '$1'" >&2; exit 2; }
  IFS=$'\t' read -r _ C_REF C_BAM BED SAMPLE _ EXTRA <<<"$row"
  read -ra EXTRA <<<"${EXTRA:-}"
  C_PRE=(-G "$C_REF" -f 0.01 -N "$SAMPLE" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z "${EXTRA[@]}"); C_POST=("$BED")
  C_DHAT=1; C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$BED")
}
cell_spec() {
  local D=scratch/g4 FA CHR1
  case "$1" in
    tiny|s1a|s8|s9|s3) slice_spec "$1" ;;
    rna-i)
      C_REF=VarDictJava/testdata/fastas/hg38.fa
      C_BAM=VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Tcov_zero_SV_DEL_candidate3.bam
      C_PRE=(-G "$C_REF" -b "$C_BAM" -R chr8:90534611-90620070 -f 0.1 -i --emit-intermediate); C_POST=()
      C_DHAT=1; C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "$C_BAM") ;;
    hg002)
      C_REF=$D/ref/giabv3/GRCh38_GIABv3.fasta; C_BAM=$D/hg002/HG002.chr20.bam
      C_PRE=(-G "$C_REF" -f 0.01 -N HG002 -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$D/hg002/eval.bed")
      C_DHAT=0; C_BIG=1; C_EST=${PROFILE_EST_MIN:-25}; C_INPUTS=("$C_REF" "$C_BAM" "$D/hg002/eval.bed") ;;
    hg008)
      C_REF=$D/ref/giabv3/GRCh38_GIABv3.fa; C_BAM="$D/hg008/HG008-T.slice.bam|$D/hg008/HG008-N.slice.bam"
      C_PRE=(-G "$C_REF" -f 0.01 -N "HG008-T|HG008-N" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$D/hg008/engine.bed")
      C_DHAT=0; C_BIG=1; C_EST=${PROFILE_EST_MIN:-20}
      C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$D/hg008/engine.bed") ;;
    wgs-scale)
      FA=$D/ref/d1vd1/GRCh38.d1.vd1.fa; C_REF=$FA; C_BAM=$D/seqc2/HCC1395.bam
      CHR1=$(awk '$1=="chr1"||$1=="1"{print $1; exit}' "$FA.fai" 2>/dev/null || true)
      C_PRE=(-G "$FA" -f 0.01 -N HCC1395 -b "$C_BAM" -R "${CHR1:-chr1}:1-10000000"); C_POST=()
      C_DHAT=0; C_BIG=1; C_EST=${PROFILE_EST_MIN:-30}; C_INPUTS=("$FA" "$FA.fai" "$C_BAM" "$C_BAM.bai") ;;
  esac
}

bin_for() {   # <cpu|heap|dhat>
  local var="PROFILE_BIN_${1^^}" base=${PROFILE_TARGET_BASE:-${CARGO_TARGET_DIR:-$ROOT_DIR/target}}
  if [ -n "${!var:-}" ]; then echo "${!var}"
  elif [ -n "${PROFILE_BIN:-}" ]; then echo "$PROFILE_BIN"
  else echo "$base/profile-$1/profiling/vardictrust"; fi
}
overridden() { local var="PROFILE_BIN_${1^^}"; [ -n "${!var:-}" ] || [ -n "${PROFILE_BIN:-}" ]; }
qcmd() { printf '%q ' "$@"; echo; }

# --- selection ---------------------------------------------------------------------------------------------
CELLS=(); for c in $ALL_CELLS; do if [ -z "$ONLY" ] || [ "$c" = "$ONLY" ]; then CELLS+=("$c"); fi; done
NEED_MODS=(cpu heap)
for c in "${CELLS[@]}"; do cell_spec "$c"; if [ "$C_DHAT" = 1 ]; then NEED_MODS+=(dhat); break; fi; done

if [ "$DRY" = 0 ]; then
  [ -x /usr/bin/time ] || { echo "profile: /usr/bin/time (GNU time) not found" >&2; exit 1; }
  for m in "${NEED_MODS[@]}"; do
    b=$(bin_for "$m")
    if ! overridden "$m"; then
      echo "profile: building profile-$m ($b)" >&2
      CARGO_TARGET_DIR="${b%/profiling/vardictrust}" cargo build --profile profiling -p vardictrust --features "profile-$m" >&2
    fi
    [ -x "$b" ] || { echo "profile: profile binary missing/not executable: $b (profile-$m)" >&2; exit 1; }
  done
  # one up-front courtesy gate (disk/RAM/maintenance window) for the longest selected cell: a refusal is final.
  maxest=1; for c in "${CELLS[@]}"; do cell_spec "$c"; [ "$C_EST" -le "$maxest" ] || maxest=$C_EST; done
  EST_MIN=$maxest bash harness/courtesy.sh true || { echo "profile: courtesy.sh refused (see above); not bypassing" >&2; exit 1; }
  GUARD=()
  if systemd-run --user --scope -p MemoryMax=14G -p MemorySwapMax=0 true >/dev/null 2>&1; then
    GUARD=(systemd-run --user --scope -p MemoryMax=14G -p MemorySwapMax=0)
  else
    echo "profile: systemd-run --user unavailable; running WITHOUT the MemoryMax=14G guard" >&2
  fi
fi

run_cell() {   # <cell>; returns 0 ok, 1 failed, 3 skipped
  local cell=$1 d="$OUT/$1" m rc=0 missing="" f mods=(cpu heap)
  cell_spec "$cell"
  for f in "${C_INPUTS[@]}"; do [ -e "$f" ] || { missing=$f; break; }; done
  if [ -n "$missing" ]; then
    if [ "$DRY" = 1 ]; then echo "NOTE $cell: input missing here: $missing"
    elif [ -n "$ONLY" ]; then echo "profile: cell '$cell': missing input $missing" >&2; return 1
    else echo "SKIP $cell: missing input $missing"; return 3; fi
  fi
  [ "$C_DHAT" = 1 ] && mods+=(dhat)
  [ "$DRY" = 1 ] && echo "PLAN $cell -> $d  (est ${C_EST} min, modalities: ${mods[*]})"
  [ "$DRY" = 1 ] || { rm -rf "${d:?}"; mkdir -p "$d"; echo "# HEAD $(git rev-parse --short HEAD 2>/dev/null || echo unknown)" > "$d/cmd.txt"; }
  for m in "${mods[@]}"; do
    local b; b=$(bin_for "$m")
    local run=(env "VARDICT_PROFILE_${m^^}=$d" "$b" "${C_PRE[@]}" -th "$TH" "${C_POST[@]}")
    local wrapped=(/usr/bin/time -v -a -o "$d/time.txt" "${run[@]}")
    if [ "$DRY" = 1 ]; then
      local note=""; if [ "$C_BIG" = 1 ] || [ "$m" = dhat ]; then note="[systemd-run MemoryMax=14G guard if available] "; fi
      printf '  [%s] EST_MIN=%s bash harness/courtesy.sh %s%s> %s\n' "$m" "$C_EST" "$note" "$(qcmd "${wrapped[@]}")" "$d/out-$m.vcf"
      continue
    fi
    if [ "$C_BIG" = 1 ] || [ "$m" = dhat ]; then wrapped=(${GUARD[@]+"${GUARD[@]}"} "${wrapped[@]}"); fi
    { echo "# modality: $m"; qcmd EST_MIN="$C_EST" bash harness/courtesy.sh "${wrapped[@]}"; } >> "$d/cmd.txt"
    echo "# modality: $m cell: $cell" >> "$d/time.txt"
    echo "profile: $cell [$m] ..." >&2
    EST_MIN="$C_EST" LC_ALL=C bash harness/courtesy.sh "${wrapped[@]}" > "$d/out-$m.vcf" || {
      rc=$?; echo "profile: FAILED $cell [$m] (exit $rc); see $d/time.txt" >&2; return 1; }
  done
  return 0
}

FAIL=0; RAN=0
for c in "${CELLS[@]}"; do
  rc=0; run_cell "$c" || rc=$?
  case $rc in 0) RAN=$((RAN + 1)) ;; 3) ;; *) FAIL=1 ;; esac
done

if [ "$DRY" = 0 ] && [ "$RAN" -gt 0 ]; then
  DIRS=(); for c in $ALL_CELLS; do
    d="$OUT/$c"; for f in cpu-collapsed.txt dhat-heap.json heap.pb.gz; do [ -e "$d/$f" ] && { DIRS+=("$d"); break; }; done
  done
  python3 harness/profile_catalog.py "${DIRS[@]}" -o "$OUT/catalog.md" >&2 || FAIL=1
  echo "profile: catalog -> $OUT/catalog.md"
fi
[ "$DRY" = 1 ] || [ "$RAN" -gt 0 ] || echo "profile: no cell ran (all inputs missing?)" >&2
exit "$FAIL"
