#!/usr/bin/env bash
set -uo pipefail
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
printf '#!/usr/bin/env bash\nprintf "%%s\\n" "$@" > "$STUB_LOG"; shift 3; exec "$@"\n' > "$T/conda"
printf '#!/usr/bin/env bash\n[ -n "${FAKE_VJ_FAIL:-}" ] && exit 7; echo fake\n' > "$T/vardict-java"
for s in teststrandbias.R var2vcf_valid.pl testsomatic.R var2vcf_paired.pl; do printf '#!/usr/bin/env bash\ncat\n' > "$T/$s"; done
chmod +x "$T"/*; export PATH="$T:$PATH" STUB_LOG="$T/args"
bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o.vcf" S -th 8 -UN -mfreq 0.2 -f0.001 >/dev/null
[ "$(grep -A1 -x -- '-th' "$T/args" | tail -1)" = 8 ] || { echo "FAIL: -th 8 not preserved"; fail=1; }
grep -qx -- '-UN' "$T/args"    || { echo "FAIL: -UN split"; fail=1; }
grep -qx -- '-mfreq' "$T/args" || { echo "FAIL: -mfreq split"; fail=1; }
grep -qx -- '0.001' "$T/args"  || { echo "FAIL: -f0.001 not split"; fail=1; }
if FAKE_VJ_FAIL=1 bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o.vcf" S >/dev/null 2>&1; then
  echo "FAIL: java crash masked (exit 0)"; fail=1; fi
[ $fail -eq 0 ] && echo "PASS test_run_oracle"; exit $fail
