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
      case "${arg:1}" in
        adaptor|chimeric|deldupvar|fisher|DP|mfreq|nmfreq|th|UN|VS)
          # Java multi-char single-dash options (CmdParser.java buildOptions) — never split.
          EXTRA_ARGS+=("$arg"); continue ;;
      esac
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

tmp_contigs=""
if [ -f "${REF}.fai" ]; then
  tmp_contigs=$(mktemp)
  awk '{print $1"\t0\t"$2}' "${REF}.fai" > "$tmp_contigs"
fi
stage_err=$(mktemp)
trap 'rm -f "$tmp_contigs" "$stage_err"' EXIT

# T3D-04 diagnostics: on engine-stage failure persist stage stderr + exit code
# + free snapshot to $OUT dir (KEEP_DIR only, one small file, no retries).
oracle_fail() {
  [ -n "${KEEP_DIR:-}" ] || return 0
  mkdir -p "$(dirname "$OUT")"
  {
    echo "oracle-stage: $1"
    echo "oracle-exit-code: $2"
    echo "--- free -h ---"
    free -h || true
    echo "--- stderr (tail 100) ---"
    tail -n 100 "$3" || true
  } > "$(dirname "$OUT")/oracle-fail-$1.log" 2>/dev/null || true
}

if [[ "$BAM" == *"|"* ]]; then  # paired/somatic
  set +e
  conda run -n vardict-oracle bash -c \
    '
    set -uo pipefail
    V2V_ARGS=()
    if [ -n "$8" ] && [ -f "$8" ]; then
      V2V_ARGS=(-G "$1" -b "$8")
    fi
    if [ "$7" = true ]; then
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:9}" "$5" | var2vcf_paired.pl -N "$3" -f "$2" "${V2V_ARGS[@]}" > "$6"
      st=("${PIPESTATUS[@]}")
      [ "${st[0]}" -eq 0 ] || { echo "oracle-stage-fail: vardict-java exit ${st[0]}" >&2; exit "${st[0]}"; }
      [ "${st[1]}" -eq 0 ] || { echo "oracle-stage-fail: var2vcf_paired.pl exit ${st[1]}" >&2; exit "${st[1]}"; }
    else
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:9}" "$5" | testsomatic.R | var2vcf_paired.pl -N "$3" -f "$2" "${V2V_ARGS[@]}" > "$6"
      st=("${PIPESTATUS[@]}")
      [ "${st[0]}" -eq 0 ] || { echo "oracle-stage-fail: vardict-java exit ${st[0]}" >&2; exit "${st[0]}"; }
      [ "${st[1]}" -eq 0 ] || { echo "oracle-stage-fail: testsomatic.R exit ${st[1]}" >&2; exit "${st[1]}"; }
      [ "${st[2]}" -eq 0 ] || { echo "oracle-stage-fail: var2vcf_paired.pl exit ${st[2]}" >&2; exit "${st[2]}"; }
    fi
    ' \
    _ "$REF" "$AF" "$SAMPLE" "$BAM" "$BED" "$OUT" "$has_fisher" "$tmp_contigs" "${EXTRA_ARGS[@]}" 2>"$stage_err"
  code=$?
  set -e
  if [ "$code" -ne 0 ]; then oracle_fail paired "$code" "$stage_err"; cat "$stage_err" >&2 || true; exit "$code"; fi
  cat "$stage_err" >&2 || true
else
  set +e
  conda run -n vardict-oracle bash -c \
    '
    set -uo pipefail
    V2V_ARGS=()
    if [ -n "$8" ] && [ -f "$8" ]; then
      V2V_ARGS=(-G "$1" -b "$8")
    fi
    if [ "$7" = true ]; then
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:9}" "$5" | var2vcf_valid.pl -N "$3" -E -f "$2" "${V2V_ARGS[@]}" > "$6"
      st=("${PIPESTATUS[@]}")
      [ "${st[0]}" -eq 0 ] || { echo "oracle-stage-fail: vardict-java exit ${st[0]}" >&2; exit "${st[0]}"; }
      [ "${st[1]}" -eq 0 ] || { echo "oracle-stage-fail: var2vcf_valid.pl exit ${st[1]}" >&2; exit "${st[1]}"; }
    else
      vardict-java -G "$1" -f "$2" -N "$3" -b "$4" -c 1 -S 2 -E 3 -g 4 "${@:9}" "$5" | teststrandbias.R | var2vcf_valid.pl -N "$3" -E -f "$2" "${V2V_ARGS[@]}" > "$6"
      st=("${PIPESTATUS[@]}")
      [ "${st[0]}" -eq 0 ] || { echo "oracle-stage-fail: vardict-java exit ${st[0]}" >&2; exit "${st[0]}"; }
      [ "${st[1]}" -eq 0 ] || { echo "oracle-stage-fail: teststrandbias.R exit ${st[1]}" >&2; exit "${st[1]}"; }
      [ "${st[2]}" -eq 0 ] || { echo "oracle-stage-fail: var2vcf_valid.pl exit ${st[2]}" >&2; exit "${st[2]}"; }
    fi
    ' \
    _ "$REF" "$AF" "$SAMPLE" "$BAM" "$BED" "$OUT" "$has_fisher" "$tmp_contigs" "${EXTRA_ARGS[@]}" 2>"$stage_err"
  code=$?
  set -e
  if [ "$code" -ne 0 ]; then oracle_fail single "$code" "$stage_err"; cat "$stage_err" >&2 || true; exit "$code"; fi
  cat "$stage_err" >&2 || true
fi
