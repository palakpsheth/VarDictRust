#!/usr/bin/env bash
set -euo pipefail
# Usage: run_cpp.sh <ref.fa> <bam | tumor.bam|normal.bam> <bed> <out.vcf> [sample] [extra VarDict args...]
# Third-arm oracle mirror of run_oracle.sh, backed by pinned vardictcpp
# (benchmarks/env/vardict-oracle.yml: vardictcpp=3=hc52dbad_0). Same argv
# contract (AF_THR env, -z translation, *|* paired form, "$@" forwarding);
# always emits VCF to OUT via -vcf on stdout (-o is numeric in cpp, not a file).
# Sharp edges (Task 1 live probes + T3B amplicon probes, never via --help/-h):
# the BED travels positionally exactly like run_oracle.sh (cpp documents a
# positional BED input; -R takes ONLY 1-based chr:start-end and -R + BED
# means the BED is silently ignored — fatal in -a mode, where cpp requires
# the primer BED). No -R expansion, no -o file (-o is numeric in cpp).
# binary picks var2vcfSingle/var2vcfPaired itself from the -b form (no perl/R).
REF=$1; BAM=$2; BED=$3; OUT=$4; SAMPLE=${5:-SAMPLE}; shift 5 || shift $#
AF=${AF_THR:-0.01}
# CPP-5 chunk pin (protocol §6/Tier-P): cpp's implicit default rewrites the
# `Seg` column, so unpinned runs are irreproducible — always pass an explicit
# --chunk (default 50000, the extract README mode-parity-table reference
# value, frozen in benchmarks/gcp-prereg.md). A caller-supplied --chunk wins,
# mirroring the -f/AF_THR rule below.
CPP_CHUNK=${CPP_CHUNK:-50000}

# Bare -z spells boolean zero-based for cpp (it takes no value: a mirrored
# "-z 1" would leak a stray positional that cpp opens as a BED file).
EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --*) EXTRA_ARGS+=("$arg") ;;
    -[a-zA-Z]*)
      case "${arg:1}" in
        adaptor|chimeric|deldupvar|fisher|DP|mfreq|nmfreq|th|UN|VS)
          EXTRA_ARGS+=("$arg"); continue ;;
      esac
      opt="${arg:0:2}"
      val="${arg:2}"
      if [ -n "$val" ]; then
        EXTRA_ARGS+=("$opt" "$val")
      elif [ "$opt" = "-z" ]; then
        EXTRA_ARGS+=("-z")
      else
        EXTRA_ARGS+=("$arg")
      fi
      ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

# Caller -f wins over AF_THR (commons-cli keeps FIRST; base -f $AF below would
# otherwise override it) — mirrors run_oracle.sh.
for ((i=0; i<${#EXTRA_ARGS[@]}-1; i++)); do
  if [ "${EXTRA_ARGS[$i]}" = "-f" ]; then
    AF="${EXTRA_ARGS[$((i+1))]}"
    break
  fi
done

# Caller --chunk wins over $CPP_CHUNK (same FIRST-wins spirit as -f above);
# otherwise append the pin so cpp never runs on its implicit default.
chunk="$CPP_CHUNK"
has_chunk=false
for ((i=0; i<${#EXTRA_ARGS[@]}; i++)); do
  case "${EXTRA_ARGS[$i]}" in
    --chunk)
      has_chunk=true
      [ $((i+1)) -lt ${#EXTRA_ARGS[@]} ] && chunk="${EXTRA_ARGS[$((i+1))]}"
      break ;;
    --chunk=*)
      has_chunk=true
      chunk="${EXTRA_ARGS[$i]#--chunk=}"
      break ;;
  esac
done
$has_chunk || EXTRA_ARGS+=(--chunk "$chunk")

# Availability gate: downstream 3-arm callers must see SKIP, not a silent 2-arm.
if ! conda run -n vardict-oracle vardictcpp --version >/dev/null 2>&1; then
  echo "SKIP: vardictcpp not available in conda env vardict-oracle" >&2
  exit 3
fi

[ -r "$BED" ] || { echo "run_cpp.sh: cannot read BED $BED" >&2; exit 1; }
# cpp requires >=4 BED columns (-c/-S/-E/-g need a gene column); pad short
# data lines with "." (the gene label is invisible in VCF output, and the
# Java oracle tolerates 3-col BEDs, so padding only restores parity of inputs).
if awk -F'\t' '!/^#/ && NF && NF<4 { missing=1; exit } END { exit !missing }' "$BED"; then
  padtmp=$(mktemp); trap 'rm -f "$padtmp"' EXIT
  awk -F'\t' 'BEGIN{OFS="\t"} /^#/ || !NF { print; next } { while (NF<4) $(NF+1)="."; print }' "$BED" > "$padtmp"
  BED="$padtmp"
fi
# Region-free BED (empty/comments only): emit a header-only VCF and exit 0,
# mirroring the oracle's empty result (compare_vcf.py + --min-variants floors
# the verdict to MISMATCH, never a vacuous MATCH).
if ! awk '!/^#/ && NF { if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/) found=1 } END { exit !found }' "$BED"; then
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$OUT"
  exit 0
fi

: > "$OUT"
# Effective chunk to the run log (compare.sh redirects this stderr to
# $OUT/cpp.err; staple_cpp_extra.sh scrapes it into extra.cpp_chunk).
echo "cpp_chunk=$chunk" >&2
# Single invocation, BED positional (java-mirror): VCF via -vcf on stdout.
conda run -n vardict-oracle vardictcpp -G "$REF" -f "$AF" -N "$SAMPLE" \
  -b "$BAM" -c 1 -S 2 -E 3 -g 4 "${EXTRA_ARGS[@]}" "$BED" -vcf > "$OUT"
