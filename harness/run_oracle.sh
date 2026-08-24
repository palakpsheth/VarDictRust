#!/usr/bin/env bash
set -euo pipefail
# Usage: run_oracle.sh <ref.fa> <bam | tumor.bam|normal.bam> <bed> <out.vcf> [sample] [extra VarDict args...]
REF=$1; BAM=$2; BED=$3; OUT=$4; SAMPLE=${5:-SAMPLE}; shift 5 || shift $#
AF=${AF_THR:-0.01}

# VarDictJava CLI doesn't support concatenated options like -f0.001; split them
EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --*) EXTRA_ARGS+=("$arg") ;;
    -[a-zA-Z]*)
      opt="${arg:0:2}"
      val="${arg:2}"
      if [ -n "$val" ]; then
        EXTRA_ARGS+=("$opt" "$val")
      elif [ "$opt" = "-z" ]; then
        # Java's -z takes an OPTIONAL numeric arg (commons-cli hasOptionalArgs,
        # CmdParser.java buildOptions). A bare trailing -z would swallow the BED
        # path as its value and crash vardict-java (silently empty VCF). Spell
        # out the default ("1") that upstream IntegrationTest DEFAULT_ARGS rely on.
        # NOTE: to disable zero-based, pass "-z0" (attached form), not "-z 0".
        EXTRA_ARGS+=("-z" "1")
      else
        EXTRA_ARGS+=("$arg")
      fi
      ;;
    *) EXTRA_ARGS+=("$arg") ;;
  esac
done

# Commons-cli keeps the FIRST value for a repeated option
# (CmdParser.java:605-608 getParsedOptionValue -> first value), so the base
# "-f $AF" below would silently override a caller-supplied case flag like
# "-f 0.001". Prefer an explicit case -f so both engines and perl see the
# same threshold (compare.sh's rust side already resolves duplicates the
# same way).
for ((i=0; i<${#EXTRA_ARGS[@]}-1; i++)); do
  if [ "${EXTRA_ARGS[$i]}" = "-f" ]; then
    AF="${EXTRA_ARGS[$((i+1))]}"
    break
  fi
done

has_fisher=false
for arg in "${EXTRA_ARGS[@]}"; do
  if [[ "$arg" == "-fisher" || "$arg" == "--fisher" ]]; then
    has_fisher=true
  fi
done

if [[ "$BAM" == *"|"* ]]; then  # paired/somatic
  conda run -n vardict-oracle bash -c \
    '
    if [ "$7" = true ]; then
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:8}" "$5" | var2vcf_paired.pl -N "$3" -f "$2" > "$6"
    else
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:8}" "$5" | testsomatic.R | var2vcf_paired.pl -N "$3" -f "$2" > "$6"
    fi
    ' \
    _ "$REF" "$AF" "$SAMPLE" "$BAM" "$BED" "$OUT" "$has_fisher" "${EXTRA_ARGS[@]}"
else
  conda run -n vardict-oracle bash -c \
    '
    if [ "$7" = true ]; then
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:8}" "$5" | var2vcf_valid.pl -N "$3" -E -f "$2" > "$6"
    else
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:8}" "$5" | teststrandbias.R | var2vcf_valid.pl -N "$3" -E -f "$2" > "$6"
    fi
    ' \
    _ "$REF" "$AF" "$SAMPLE" "$BAM" "$BED" "$OUT" "$has_fisher" "${EXTRA_ARGS[@]}"
fi
