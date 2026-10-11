#!/usr/bin/env bash
set -euo pipefail
# Usage: run_cpp.sh <ref.fa> <bam | tumor.bam|normal.bam> <bed> <out.vcf> [sample] [extra VarDict args...]
# Third-arm oracle mirror of run_oracle.sh, backed by pinned vardictcpp
# (benchmarks/env/vardict-oracle.yml: vardictcpp=3=hc52dbad_0). Same argv
# contract (AF_THR env, -z translation, *|* paired form, "$@" forwarding);
# always emits VCF to OUT via -vcf on stdout (-o is numeric in cpp, not a file).
# Sharp edges (Task 1 live probes, never re-probe via --help/-h): -R takes ONLY
# 1-based chr:start-end (BED paths rejected), so the BED is expanded to region
# strings (start+1) with one invocation per region (last -R wins); the single
# binary picks var2vcfSingle/var2vcfPaired itself from the -b form (no perl/R).
REF=$1; BAM=$2; BED=$3; OUT=$4; SAMPLE=${5:-SAMPLE}; shift 5 || shift $#
AF=${AF_THR:-0.01}

# Arg-shifting mirrors run_oracle.sh (VarDictJava CLI has no concatenated
# -f0.001; multi-char single-dash options are never split; bare -z spells 1).
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
        EXTRA_ARGS+=("-z" "1")
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

# Availability gate: downstream 3-arm callers must see SKIP, not a silent 2-arm.
if ! conda run -n vardict-oracle vardictcpp --version >/dev/null 2>&1; then
  echo "SKIP: vardictcpp not available in conda env vardict-oracle" >&2
  exit 3
fi

# BED (0-based, half-open) -> 1-based inclusive region strings for -R.
REGIONS=()
while IFS=$'\t' read -r chrom s e _ || [ -n "${chrom:-}" ]; do
  case "${chrom:-}" in ''|\#*) continue;; esac
  [[ "$s" =~ ^[0-9]+$ && "$e" =~ ^[0-9]+$ ]] || continue
  REGIONS+=("$chrom:$((s+1))-$e")
done < "$BED"

: > "$OUT"
if [ "${#REGIONS[@]}" -gt 0 ]; then
  tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
  for R in "${REGIONS[@]}"; do
    conda run -n vardict-oracle vardictcpp -G "$REF" -f "$AF" -N "$SAMPLE" \
      -b "$BAM" -c 1 -S 2 -E 3 -g 4 "${EXTRA_ARGS[@]}" -R "$R" -vcf > "$tmp"
    if [ -s "$OUT" ]; then
      awk '!/^#/' "$tmp" >> "$OUT"  # header already written; awk exits 0 on empty
    else
      cat "$tmp" > "$OUT"
    fi
  done
fi
