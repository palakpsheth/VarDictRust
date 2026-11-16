#!/usr/bin/env bash
# Usage: concordance_pair.sh <name> <java.vcf> <rust.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor> [cpp.vcf]
# Gate (D-27): identical 'None'-row counts AND identical tp/tp-baseline/fp/fn site lists. Not a dF1 tolerance.
# Optional 9th arg adds the advisory vardictcpp third arm (pinned vardictcpp 3
# hc52dbad_0 bioconda): a second concordance.sh vcfeval run with identical
# truth/BED/sdf. cpp F1 lands in extra.json (cpp_f1 + cpp_version) and stdout
# as a rank line; the gate and exit code stay java==rust only (Task-3 ruling:
# header/END-key deltas are known-normalized, cpp rank never gated).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/portable.sh"
source "$SCRIPT_DIR/lib/concordance_cpp.sh"

if ! command -v rtg &>/dev/null || ! command -v bcftools &>/dev/null; then
  export PATH="$G4_TOOLS_BIN:$PATH"
fi

NAME=$1; J=$2; R=$3; TRUTH=$4; BED=$5; SDF=$6; OUT=$7; TS=$8; CPP_VCF=${9:-}
bash harness/concordance.sh "$J" "$TRUTH" "$BED" "$SDF" "$OUT/java" "$TS" >/dev/null
bash harness/concordance.sh "$R" "$TRUTH" "$BED" "$SDF" "$OUT/rust" "$TS" >/dev/null

run_pair_cpp_arm() {
  local arm="$1" vcf_in="$2"
  bash harness/concordance.sh "$vcf_in" "$TRUTH" "$BED" "$SDF" "$OUT/$arm" "$TS" >"$OUT/${arm}.stdout"
}

# Third arm: vardictcpp (advisory only via lib/concordance_cpp.sh; sets CPP_F1 / CPP_F1_NORM / CPP_NOTE)
concordance_cpp_run_third_arm "$CPP_VCF" "$TRUTH" "$J" "$OUT" run_pair_cpp_arm

res=EQUAL; cmp -s "$OUT/java/summary.json" "$OUT/rust/summary.json" || res=DIFFERENT
for f in tp tp-baseline fp fn; do cmp -s "$OUT/java/$f.sites" "$OUT/rust/$f.sites" || res=DIFFERENT; done
# Emits extra.json (java/rust summaries + cpp_f1 / cpp_f1_normalized / cpp_version or cpp_skip_reason)
concordance_cpp_write_extra "$OUT" "$CPP_F1" "$CPP_NOTE" "$CPP_F1_NORM"
python3 harness/record_evidence.py --input "$TRUTH" --input "$BED" "$NAME" "$res" "$OUT" -- "$0" "$@"
[ "$res" = EQUAL ]
