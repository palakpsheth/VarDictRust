#!/usr/bin/env bash
set -uo pipefail
# Test for harness/run_cpp.sh (vardictcpp third-arm oracle mirror).
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0

# 1. Missing env -> SKIP exit 3, never a silent 2-arm (stub conda, subshell PATH).
# Runs unconditionally even on hosts without the live vardict-oracle env.
S="$T/stub"; mkdir -p "$S"
printf '#!/usr/bin/env bash\necho "stub conda: no env" >&2; exit 1\n' > "$S/conda"
chmod +x "$S/conda"
skip_out=$(PATH="$S:$PATH" bash harness/run_cpp.sh r.fa b.bam x.bed "$T/o.vcf" 2>&1); rc=$?
[ "$rc" -eq 3 ] || { echo "FAIL: missing env exit=$rc, want 3"; fail=1; }
[[ "$skip_out" == SKIP:*vardictcpp* ]] || { echo "FAIL: missing env msg: $skip_out"; fail=1; }

# 2. Availability gate: without the live conda env, skip remaining live tests cleanly (exit 0 if stub test passed).
if ! conda run -n vardict-oracle vardictcpp --version >/dev/null 2>&1; then
  echo "SKIP: vardictcpp not available in conda env vardict-oracle"
  [ $fail -eq 0 ] && echo "PASS test_run_cpp"
  exit $fail
fi

# 3. Version pin (Task 1: vardictcpp=3=hc52dbad_0).
ver=$(conda run -n vardict-oracle vardictcpp --version 2>&1)
[[ "$ver" == *"vardictcpp 3"* ]] || { echo "FAIL: --version pin mismatch: $ver"; fail=1; }
conda list -n vardict-oracle 2>/dev/null | grep -q 'vardictcpp.*hc52dbad_0' \
  || { echo "FAIL: conda pin vardictcpp=3=hc52dbad_0 missing"; fail=1; }

# 4. Micro fixture -> non-empty VCF (bcftools record view, never head -1).
bash harness/run_cpp.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam \
  fixtures/micro/micro.bed "$T/cpp.vcf" >/dev/null 2>"$T/cpp.err" \
  || { echo "FAIL: run_cpp.sh exit $?"; fail=1; }
conda run -n vardict-oracle bcftools view -H "$T/cpp.vcf" > "$T/h.txt" 2>/dev/null \
  || { echo "FAIL: bcftools cannot parse cpp.vcf"; fail=1; }
grep -qm1 . "$T/h.txt" || { echo "FAIL: cpp.vcf has no records"; fail=1; }

# 5. <4-col BED padding (run_cpp.sh pads 3-col BED with synthetically added 4th column ".").
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' fixtures/micro/micro.bed > "$T/micro_3col.bed"
bash harness/run_cpp.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam \
  "$T/micro_3col.bed" "$T/cpp_3col.vcf" >/dev/null 2>"$T/cpp_3col.err" \
  || { echo "FAIL: run_cpp.sh 3-col BED exit $?"; fail=1; }
conda run -n vardict-oracle bcftools view -H "$T/cpp_3col.vcf" > "$T/h3.txt" 2>/dev/null \
  || { echo "FAIL: bcftools cannot parse cpp_3col.vcf"; fail=1; }
cmp -s "$T/h.txt" "$T/h3.txt" || { echo "FAIL: 3-col BED records differ from 4-col BED records"; fail=1; }

[ $fail -eq 0 ] && echo "PASS test_run_cpp"; exit $fail
