#!/usr/bin/env bash
# Usage: profile.sh baseline|optimized [--only <cell>] [--dry-run]
#   Profile VarDictRust (Rust engine only, no oracle) on a fixed set of named cells with the three profiler builds
#   of crates/vardictrust (features profile-cpu / profile-heap / profile-dhat) and rank the result:
#     scratch/profiles/<tag>/<cell>/{cpu-*,heap*,dhat-*,time.txt,cmd.txt,out-*}  +  scratch/profiles/<tag>/catalog.md
#   (catalog = harness/profile_catalog.py; SVG flamegraphs per cell are linked from it.)
#
# Cells (commands copied from existing harness files, never invented):
#   tiny s1a s1b        germline slices, harness/slices.tsv rows (run_slice.sh columns)   [dhat]
#   s8 s9               indel/SV slices, harness/slices.tsv                               [dhat]
#   s3 s4 s5            somatic tumor|normal pairs, harness/slices.tsv                    [dhat]
#   s6                  single-amplicon 8-col BED (-a 10:0.95), harness/slices.tsv        [dhat]
#   s7                  pileup + amplicon (-a 10:0.95 -p), harness/slices.tsv             [dhat]
#   amplicon            overlapping 2-amplicon 8-col BED (amplicon-tier2.json cmd)        [dhat]
#   crispr              extended m4: -R chr20:50454281-50455444 -J 50454941 -j 25         [dhat]
#   rna-i               extended m1: -R chr8:90534611-90620070 -f 0.1 -i                  [dhat]
#   hg005               multimodal HG005_1A: 300x Han Chinese WGS slice (-f 0.05)         [dhat]
#   seqc2-fisher        multimodal SEQC2_2A: somatic WES slice with --fisher              [dhat]
#   seqc2-slices        SEQC2 HCC1395|BL somatic slices (slices.engine.bed)               [no dhat]
#   ext-*               extended m5..m16 CLI flag modes (run_extended_coverage.sh)        [dhat]
#   hg002               HG002 chr20:10M-30M slice (run_hg002.sh arguments)                [no dhat]
#   hg008               HG008 T|N slice, engine.bed (benchmarks/hg008-v04-guide.md cmd)   [no dhat]
#   hg008-cram          HG008 T|N CRAM twin, engine.bed (cram_twin.sh, RUST_MIN_STACK=64M)[no dhat]
#   wgs-scale           scratch/g4/seqc2/HCC1395.bam, bounded 10 Mb -R window at -th 8    [no dhat]
#   hg002-sub           opt-in (--only): first 83 regions (<10.6M, 594.5 kb) of hg002/eval.bed  [dhat only]
#   hg008-sub           opt-in (--only): first 60 regions (774.7 kb) of hg008/engine.bed        [dhat only]
#   hg008-cram-sub      opt-in (--only): first 60 regions (774.7 kb) of hg008/engine.bed CRAM   [dhat only]
#   seqc2-slices-sub    opt-in (--only): first 60 regions of seqc2/slices.engine.bed            [dhat only]
#   seqc2-exome         opt-in (--only): full SEQC2 HCC1395|BL exome (exome.engine.bed)         [no dhat]
#   cells whose inputs are missing are SKIPped (an error if named via --only). dhat only ever runs on slice cells
#   and the opt-in *-sub cells (its live-block map costs GBs: the systemd-oomd failure mode); HG002/HG008/SEQC2/WGS get
#   cpu + heap only. Set PROFILE_LOOPS=N (default 1) to repeat the cpu modality N times on slice cells (C_BIG=0) and merge collapsed stacks.
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

DEFAULT_CELLS="tiny s1a s1b s8 s9 s3 s4 s5 s6 s7 amplicon crispr rna-i hg005 seqc2-fisher seqc2-slices ext-unique ext-adaptor ext-norealign ext-nosv ext-nbase ext-dedup ext-shift3 ext-mapq ext-extend ext-segmap ext-minmatch ext-deldupvar hg002 hg008 hg008-cram wgs-scale"
ALL_CELLS="$DEFAULT_CELLS hg002-sub hg008-sub hg008-cram-sub seqc2-slices-sub seqc2-exome"
if [ -n "$ONLY" ]; then
  case " $ALL_CELLS " in *" $ONLY "*) ;; *) echo "profile: unknown cell '$ONLY' (known: $ALL_CELLS)" >&2; exit 2 ;; esac
fi

TH=${PROFILE_THREADS:-8}
[[ "$TH" =~ ^[0-9]+$ ]] && [ "$TH" -ge 1 ] && [ "$TH" -le 8 ] \
  || { echo "profile: PROFILE_THREADS=$TH out of range 1..8 (Host Policy A3 cap)" >&2; exit 2; }
LOOPS=${PROFILE_LOOPS:-1}
[[ "$LOOPS" =~ ^[0-9]+$ ]] && [ "$LOOPS" -ge 1 ] \
  || { echo "profile: PROFILE_LOOPS=$LOOPS must be >= 1" >&2; exit 2; }

OUT_ROOT=${PROFILE_OUT:-scratch/profiles}
OUT="$OUT_ROOT/$TAG"
EVIDENCE=$(realpath -m benchmarks/evidence)
case "$(realpath -m "$OUT")/" in
  "$EVIDENCE"/*) echo "profile: refusing to write under benchmarks/evidence/ ($OUT)" >&2; exit 2 ;;
esac

# --- cells -----------------------------------------------------------------------------------------------
# cell_spec <name> sets C_REF, C_BAM (a|b = tumor|normal), C_PRE (args before -th), C_POST (BED positional),
# C_BIG (0/1), C_EST (courtesy EST_MIN minutes), C_INPUTS (files that must exist), C_MODS, C_ENV, C_GEN_BED.
slice_spec() {   # same columns run_slice.sh reads from harness/slices.tsv; compare.sh's Rust argv shape
  local row EXTRA SAMPLE BED; row=$(awk -F'\t' -v id="$1" '$1==id' harness/slices.tsv)
  [ -n "$row" ] || { echo "profile: no harness/slices.tsv row '$1'" >&2; exit 2; }
  IFS=$'\t' read -r _ C_REF C_BAM BED SAMPLE _ EXTRA <<<"$row"
  read -ra EXTRA <<<"${EXTRA:-}"
  C_PRE=(-G "$C_REF" -f 0.01 -N "$SAMPLE" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z "${EXTRA[@]}"); C_POST=("$BED")
  C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$BED"); C_MODS=(cpu heap dhat)
}
ext_slice_spec() {   # <base_slice_id> <extra_flags...>; uses run_extended_coverage.sh default -f 0.01 (ignores slices.tsv EXTRA)
  local base=$1 row SAMPLE BED; shift
  row=$(awk -F'\t' -v id="$base" '$1==id' harness/slices.tsv)
  [ -n "$row" ] || { echo "profile: no harness/slices.tsv row '$base'" >&2; exit 2; }
  IFS=$'\t' read -r _ C_REF C_BAM BED SAMPLE _ _ <<<"$row"
  C_PRE=(-G "$C_REF" -f 0.01 -N "$SAMPLE" -b "$C_BAM" "$@" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$BED")
  C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$BED"); C_MODS=(cpu heap dhat)
}
hg008_spec() {   # <cell> <bam|cram> <full|sub>
  local cell=$1 ext=$2 mode=$3 D=scratch/g4
  C_REF=$D/ref/giabv3/GRCh38_GIABv3.fa; C_BAM="$D/hg008/HG008-T.slice.$ext|$D/hg008/HG008-N.slice.$ext"
  [ "$ext" = cram ] && C_ENV=(RUST_MIN_STACK=67108864)
  C_PRE=(-G "$C_REF" -f 0.01 -N "HG008-T|HG008-N" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z)
  C_BIG=1
  if [ "$mode" = sub ]; then
    C_SUB_SRC="$D/hg008/engine.bed"; C_SUB_CMD='head -n 60'; C_POST=("$OUT/$cell/sub.bed")
    C_EST=${PROFILE_EST_MIN:-30}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$C_SUB_SRC"); C_MODS=(dhat)
  else
    C_POST=("$D/hg008/engine.bed"); C_EST=${PROFILE_EST_MIN:-20}
    C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$D/hg008/engine.bed"); C_MODS=(cpu heap)
  fi
}
cell_spec() {
  local D=scratch/g4 FA CHR1
  C_MODS=(cpu heap); C_SUB_SRC=""; C_SUB_CMD=""; C_GEN_BED=""; C_ENV=()
  case "$1" in
    tiny|s1a|s1b|s8|s9|s3|s4|s5|s6|s7) slice_spec "$1" ;;
    amplicon)
      C_REF=$D/amplicon/ref.fa; C_BAM=$D/amplicon/amplicon.bam
      C_PRE=(-G "$C_REF" -f 0.01 -N SAMPLE -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z -a 10:0.95); C_POST=("$D/amplicon/correct_8columns_amplicon.bed")
      C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "$C_BAM" "$D/amplicon/correct_8columns_amplicon.bed"); C_MODS=(cpu heap dhat) ;;
    crispr)
      C_REF=VarDictJava/testdata/fastas/hg19.fa
      C_BAM=VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CRISPR.bam
      C_PRE=(-G "$C_REF" -b "$C_BAM" -R chr20:50454281-50455444 -f 0.001 -J 50454941 -j 25 -c 1 -S 2 -E 3 -g 4 --emit-intermediate); C_POST=()
      C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "$C_BAM"); C_MODS=(cpu heap dhat) ;;
    rna-i)
      C_REF=VarDictJava/testdata/fastas/hg38.fa
      C_BAM=VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Tcov_zero_SV_DEL_candidate3.bam
      C_PRE=(-G "$C_REF" -b "$C_BAM" -R chr8:90534611-90620070 -f 0.1 -i --emit-intermediate); C_POST=()
      C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "$C_BAM"); C_MODS=(cpu heap dhat) ;;
    hg005)
      C_REF=$D/ref/giabv3/GRCh38_GIABv3.fasta; C_BAM=$D/multimodal/hg005_s4a.bam
      C_GEN_BED=$'chr3\t40000000\t40010000\n'
      C_PRE=(-G "$C_REF" -f 0.05 -N HG005 -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$OUT/$1/sub.bed")
      C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "$C_BAM"); C_MODS=(cpu heap dhat) ;;
    seqc2-fisher)
      C_REF=$D/ref/d1vd1/GRCh38.d1.vd1.fa
      C_BAM="$D/multimodal/seqc2_s3a_tumor.bam|$D/multimodal/seqc2_s3a_normal.bam"
      C_GEN_BED=$'chr7\t140400000\t140420000\n'
      C_PRE=(-G "$C_REF" -f 0.05 -N "TUMOR|NORMAL" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z --fisher); C_POST=("$OUT/$1/sub.bed")
      C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}"); C_MODS=(cpu heap dhat) ;;
    seqc2-slices)
      C_REF=$D/ref/d1vd1/GRCh38.d1.vd1.fa; C_BAM="$D/seqc2/HCC1395.bam|$D/seqc2/HCC1395BL.bam"
      C_PRE=(-G "$C_REF" -f 0.01 -N "HCC1395|HCC1395BL" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$D/seqc2/slices.engine.bed")
      C_BIG=1; C_EST=${PROFILE_EST_MIN:-10}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$D/seqc2/slices.engine.bed") ;;
    seqc2-slices-sub)
      C_REF=$D/ref/d1vd1/GRCh38.d1.vd1.fa; C_BAM="$D/seqc2/HCC1395.bam|$D/seqc2/HCC1395BL.bam"
      C_SUB_SRC="$D/seqc2/slices.engine.bed"; C_SUB_CMD='head -n 60'
      C_PRE=(-G "$C_REF" -f 0.01 -N "HCC1395|HCC1395BL" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$OUT/$1/sub.bed")
      C_BIG=1; C_EST=${PROFILE_EST_MIN:-15}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$C_SUB_SRC"); C_MODS=(dhat) ;;
    seqc2-exome)
      C_REF=$D/ref/d1vd1/GRCh38.d1.vd1.fa; C_BAM="$D/seqc2/HCC1395.bam|$D/seqc2/HCC1395BL.bam"
      C_PRE=(-G "$C_REF" -f 0.01 -N "HCC1395|HCC1395BL" -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$D/seqc2/exome.engine.bed")
      C_BIG=1; C_EST=${PROFILE_EST_MIN:-30}; C_INPUTS=("$C_REF" "${C_BAM%%|*}" "${C_BAM##*|}" "$D/seqc2/exome.engine.bed") ;;
    ext-unique)    ext_slice_spec s1a -u -UN ;;
    ext-adaptor)   ext_slice_spec s1a --adaptor AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ;;
    ext-norealign) ext_slice_spec s1a -k 0 ;;
    ext-nosv)      ext_slice_spec s1a -U ;;
    ext-nbase)     ext_slice_spec s1a -K ;;
    ext-dedup)     ext_slice_spec tiny -t ;;
    ext-shift3)    ext_slice_spec tiny -3 ;;
    ext-mapq)      ext_slice_spec tiny -O 55 ;;
    ext-extend)    ext_slice_spec tiny -x 50 ;;
    ext-segmap)
      C_REF=VarDictJava/testdata/fastas/GRCh37.fa
      C_BAM=VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_20156208_20234556.bam
      C_GEN_BED=$'chr20\t20178900\t20179200\ttestbed\t0\t+\t20178900\t20179200\t0\t20178900\t20179200\ttestbed\n'
      C_PRE=(-G "$C_REF" -f 0.01 -N SAMPLE -b "$C_BAM" -c 1 -S 2 -E 3 -s 10 -e 11 -g 4 -z); C_POST=("$OUT/$1/sub.bed")
      C_BIG=0; C_EST=${PROFILE_EST_MIN:-5}; C_INPUTS=("$C_REF" "$C_BAM"); C_MODS=(cpu heap dhat) ;;
    ext-minmatch)  ext_slice_spec tiny -M 50 ;;
    ext-deldupvar) ext_slice_spec tiny --deldupvar ;;
    hg002)
      C_REF=$D/ref/giabv3/GRCh38_GIABv3.fasta; C_BAM=$D/hg002/HG002.chr20.bam
      C_PRE=(-G "$C_REF" -f 0.01 -N HG002 -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$D/hg002/eval.bed")
      C_BIG=1; C_EST=${PROFILE_EST_MIN:-25}; C_INPUTS=("$C_REF" "$C_BAM" "$D/hg002/eval.bed") ;;
    hg002-sub)
      C_REF=$D/ref/giabv3/GRCh38_GIABv3.fasta; C_BAM=$D/hg002/HG002.chr20.bam
      C_SUB_SRC="$D/hg002/eval.bed"; C_SUB_CMD='awk "\$2<10600000"'
      C_PRE=(-G "$C_REF" -f 0.01 -N HG002 -b "$C_BAM" -c 1 -S 2 -E 3 -g 4 -z); C_POST=("$OUT/$1/sub.bed")
      C_BIG=1; C_EST=${PROFILE_EST_MIN:-15}; C_INPUTS=("$C_REF" "$C_BAM" "$C_SUB_SRC"); C_MODS=(dhat) ;;
    hg008)          hg008_spec "$1" bam  full ;;
    hg008-sub)      hg008_spec "$1" bam  sub  ;;
    hg008-cram)     hg008_spec "$1" cram full ;;
    hg008-cram-sub) hg008_spec "$1" cram sub  ;;
    wgs-scale)
      FA=$D/ref/d1vd1/GRCh38.d1.vd1.fa; C_REF=$FA; C_BAM=$D/seqc2/HCC1395.bam
      CHR1=$(awk '$1=="chr1"||$1=="1"{print $1; exit}' "$FA.fai" 2>/dev/null || true)
      C_PRE=(-G "$FA" -f 0.01 -N HCC1395 -b "$C_BAM" -R "${CHR1:-chr1}:1-10000000"); C_POST=()
      C_BIG=1; C_EST=${PROFILE_EST_MIN:-30}; C_INPUTS=("$FA" "$FA.fai" "$C_BAM" "$C_BAM.bai") ;;
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
CELLS=()
if [ -n "$ONLY" ]; then CELLS=("$ONLY")
else for c in $DEFAULT_CELLS; do CELLS+=("$c"); done; fi
NEED_MODS=()
for m in cpu heap dhat; do
  for c in "${CELLS[@]}"; do
    cell_spec "$c"
    case " ${C_MODS[*]} " in *" $m "*) NEED_MODS+=("$m"); break ;; esac
  done
done

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
  local cell=$1 d="$OUT/$1" m rc=0 missing="" f mods=() loops
  cell_spec "$cell"; mods=("${C_MODS[@]}")
  loops=$LOOPS; [ "$C_BIG" = 1 ] && loops=1
  for f in "${C_INPUTS[@]}"; do [ -e "$f" ] || { missing=$f; break; }; done
  if [ -n "$missing" ]; then
    if [ "$DRY" = 1 ]; then echo "NOTE $cell: input missing here: $missing"
    elif [ -n "$ONLY" ]; then echo "profile: cell '$cell': missing input $missing" >&2; return 1
    else echo "SKIP $cell: missing input $missing"; return 3; fi
  fi
  [ "$DRY" = 1 ] && echo "PLAN $cell -> $d  (est ${C_EST} min, modalities: ${mods[*]}, loops: $loops)"
  if [ "$DRY" = 0 ]; then
    rm -rf "${d:?}"; mkdir -p "$d"
    echo "# HEAD $(git rev-parse --short HEAD 2>/dev/null || echo unknown)" > "$d/cmd.txt"
    if [ -n "$C_SUB_SRC" ]; then eval "$C_SUB_CMD" < "$C_SUB_SRC" > "$d/sub.bed"; fi
    if [ -n "$C_GEN_BED" ]; then printf '%s' "$C_GEN_BED" > "$d/sub.bed"; fi
  fi
  for m in "${mods[@]}"; do
    local b; b=$(bin_for "$m")
    local run=(env ${C_ENV[@]+"${C_ENV[@]}"} "VARDICT_PROFILE_${m^^}=$d" "$b" "${C_PRE[@]}" -th "$TH" "${C_POST[@]}")
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
    if [ "$m" = cpu ] && [ "$loops" -gt 1 ]; then
      for ((i=1; i<=loops; i++)); do
        mkdir -p "$d/iter-$i"
        local irun=(env ${C_ENV[@]+"${C_ENV[@]}"} "VARDICT_PROFILE_CPU=$d/iter-$i" "$b" "${C_PRE[@]}" -th "$TH" "${C_POST[@]}")
        local iwrap=(/usr/bin/time -v -a -o "$d/time.txt" "${irun[@]}")
        if [ "$C_BIG" = 1 ]; then iwrap=(${GUARD[@]+"${GUARD[@]}"} "${iwrap[@]}"); fi
        EST_MIN="$C_EST" LC_ALL=C bash harness/courtesy.sh "${iwrap[@]}" > "$d/out-$m.vcf" || {
          rc=$?; echo "profile: FAILED $cell [$m iter $i] (exit $rc); see $d/time.txt" >&2; return 1; }
      done
      awk '{n=$NF; $NF=""; sub(/[[:space:]]+$/,""); c[$0]+=n} END{for(k in c) print k, c[k]}' "$d"/iter-*/cpu-collapsed.txt | LC_ALL=C sort > "$d/cpu-collapsed.txt"
      python3 -c 'import glob,json,sys; d=sys.argv[1]; ms=[json.load(open(p)) for p in glob.glob(f"{d}/iter-*/cpu-meta.json")]; json.dump({"ticks_total":sum(m.get("ticks_total",0) for m in ms),"samples_kept":sum(m.get("samples_kept",0) for m in ms),"loops":len(ms)},open(f"{d}/cpu-meta.json","w"))' "$d"
      for f in cpu-profile.pb cpu-flamegraph.svg; do [ -e "$d/iter-1/$f" ] && cp -f "$d/iter-1/$f" "$d/$f" || true; done
    else
      EST_MIN="$C_EST" LC_ALL=C bash harness/courtesy.sh "${wrapped[@]}" > "$d/out-$m.vcf" || {
        rc=$?; echo "profile: FAILED $cell [$m] (exit $rc); see $d/time.txt" >&2; return 1; }
    fi
  done
  return 0
}

FAIL=0; RAN=0
for c in "${CELLS[@]}"; do
  rc=0; run_cell "$c" || rc=$?
  case $rc in 0) RAN=$((RAN + 1)) ;; 3) ;; *) FAIL=1; break ;; esac
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
