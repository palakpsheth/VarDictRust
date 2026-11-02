#!/usr/bin/env bash
set -uo pipefail
# Test for harness/run_oracle.sh and harness/compare.sh platform guard (Issue B).
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
export GOLDEN_STORE_DIR="$T/store"

mkdir -p "$T/bin"
cat > "$T/bin/uname" <<'EOF'
#!/usr/bin/env bash
case "${1:-}" in
  -s) echo "${FAKE_OS:-Darwin}" ;;
  -m) echo "${FAKE_ARCH:-arm64}" ;;
  *) echo "${FAKE_OS:-Darwin} ${FAKE_ARCH:-arm64}" ;;
esac
EOF
printf '#!/usr/bin/env bash\nprintf "%%s\\n" "$@" > "$STUB_LOG"; shift 3; exec "$@"\n' > "$T/bin/conda"
printf '#!/usr/bin/env bash\necho fake\n' > "$T/bin/vardict-java"
for s in teststrandbias.R var2vcf_valid.pl testsomatic.R var2vcf_paired.pl; do
  printf '#!/usr/bin/env bash\ncat\n' > "$T/bin/$s"
done
cat > "$T/bin/docker" <<'EOF'
#!/usr/bin/env bash
printf "%s\n" "$@" > "$DOCKER_LOG"
exit 0
EOF
chmod +x "$T/bin"/*

export STUB_LOG="$T/args" DOCKER_LOG="$T/docker_args"

# 1. Darwin/arm64 without ORACLE_CONTAINER must fail fast (exit 2) with ORACLE-PLATFORM on stderr.
err1=$( env -u ORACLE_CONTAINER PATH="$T/bin:$PATH" bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o1.vcf" S 2>&1 )
rc1=$?
[ "$rc1" -eq 2 ] || { echo "FAIL1: expected exit 2 on Darwin/arm64, got $rc1"; fail=1; }
grep -q 'ORACLE-PLATFORM: unsupported oracle host Darwin/arm64 (R long double divergence); set ORACLE_CONTAINER=1 to run the linux/amd64 oracle container.' <<<"$err1" \
  || { echo "FAIL1: missing ORACLE-PLATFORM message on stderr: $err1"; fail=1; }

# 2. Linux/aarch64 without ORACLE_CONTAINER must fail fast (exit 2) with ORACLE-PLATFORM on stderr.
err2=$( env -u ORACLE_CONTAINER FAKE_OS=Linux FAKE_ARCH=aarch64 PATH="$T/bin:$PATH" bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o2.vcf" S 2>&1 )
rc2=$?
[ "$rc2" -eq 2 ] || { echo "FAIL2: expected exit 2 on Linux/aarch64, got $rc2"; fail=1; }
grep -q 'ORACLE-PLATFORM: unsupported oracle host Linux/aarch64 (R long double divergence); set ORACLE_CONTAINER=1 to run the linux/amd64 oracle container.' <<<"$err2" \
  || { echo "FAIL2: missing ORACLE-PLATFORM message on stderr: $err2"; fail=1; }

# 3. compare.sh on Darwin/arm64 without ORACLE_CONTAINER must also fail fast (exit 2) with ORACLE-PLATFORM.
err3=$( env -u ORACLE_CONTAINER PATH="$T/bin:$PATH" bash harness/compare.sh r.fa b.bam x.bed S 2>&1 )
rc3=$?
[ "$rc3" -eq 2 ] || { echo "FAIL3: expected compare.sh exit 2 on Darwin/arm64, got $rc3"; fail=1; }
grep -q 'ORACLE-PLATFORM: unsupported oracle host Darwin/arm64 (R long double divergence); set ORACLE_CONTAINER=1 to run the linux/amd64 oracle container.' <<<"$err3" \
  || { echo "FAIL3: compare.sh missing ORACLE-PLATFORM message: $err3"; fail=1; }

# 4. Darwin/arm64 with ORACLE_CONTAINER=1 dispatches via docker --platform linux/amd64.
ORACLE_CONTAINER=1 PATH="$T/bin:$PATH" bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o4.vcf" S -th 4 >/dev/null 2>&1
rc4=$?
[ "$rc4" -eq 0 ] || { echo "FAIL4: expected exit 0 with ORACLE_CONTAINER=1, got $rc4"; fail=1; }
[ "$(grep -A1 -x -- '--platform' "$DOCKER_LOG" | tail -1)" = "linux/amd64" ] \
  || { echo "FAIL4: docker not invoked with --platform linux/amd64"; fail=1; }

# 5. Linux/x86_64 succeeds directly without ORACLE-PLATFORM error.
err5=$( env -u ORACLE_CONTAINER FAKE_OS=Linux FAKE_ARCH=x86_64 PATH="$T/bin:$PATH" bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o5.vcf" S -th 4 2>&1 )
rc5=$?
[ "$rc5" -eq 0 ] || { echo "FAIL5: expected exit 0 on Linux/x86_64, got $rc5 ($err5)"; fail=1; }
! grep -q 'ORACLE-PLATFORM' <<<"$err5" || { echo "FAIL5: unexpected ORACLE-PLATFORM on Linux/x86_64"; fail=1; }

[ "$fail" -eq 0 ] && echo "PASS test_oracle_platform"
exit "$fail"
