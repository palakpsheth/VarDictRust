#!/usr/bin/env bash
# Nightly heavy-evidence regen for the VarDictRust paper study.
# Scheduling vehicle is ONE-SHOT ONLY (no timer, no re-arm):
#   sudo systemd-run --on-calendar="2026-10-16 22:00" \
#     --unit=vardict-regen-nightly --collect \
#     --description="VarDict heavy evidence regen" \
#     /bin/bash /mnt/media_backups/VarDictRust/docs/validation/provenance/regen-nightly-runbook.sh
# The hardcoded `[ "$(date +%F)" = "2026-10-16" ]` guard below is intentional:
# it makes any accidental re-fire (stale timer, manual re-run) abort instead of
# re-stamping evidence. To run on another date, edit the guard deliberately.
# Runs in .worktrees/eng-v102 (branch perf/bounded-amplicon-memory).
# NEVER hand-edits evidence: every JSON is stamped by harness/record_evidence.py
# (or a self-stamping driver) only after asserting MATCH/EQUAL/PASS from live output.
# Host restore (oomd + media stack) runs in a trap on EVERY exit path.
set -euo pipefail
WT=/mnt/media_backups/VarDictRust/.worktrees/eng-v102
LOG=/mnt/media_backups/VarDictRust/scratch/regen-nightly.log
LOCK=/mnt/media_backups/VarDictRust/scratch/.regen-nightly.lock  # palak-owned dir: never root-stranded like /tmp was on 2026-10-16
VENV_PY=/home/palak/Documents/venv/bin/python
LIFECYCLE=/home/palak/Documents/media_stack_lifecycle.py
OOMD_WAS_ACTIVE="unknown"
STACK_DOWN=0
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
bring_up() {
  rc=$?
  set +e
  log "TRAP: restoring host (oomd + media stack)"
  if [ "$OOMD_WAS_ACTIVE" = "active" ]; then sudo systemctl start systemd-oomd 2>&1 | tee -a "$LOG"; fi
  if [ "$STACK_DOWN" = 1 ]; then
    "$VENV_PY" "$LIFECYCLE" up 2>&1 | tee -a "$LOG"
    "$VENV_PY" "$LIFECYCLE" status 2>&1 | tail -n 3 | tee -a "$LOG"
  fi
  log "TRAP DONE (script exit was $rc)"
}
[ "$(id -un)" = "palak" ] || { echo "FATAL: must run as palak"; exit 8; }
exec 9>"$LOCK"; flock -n 9 || { echo "another regen run holds $LOCK"; exit 7; }
trap bring_up EXIT
export PATH="/home/palak/miniconda3/envs/g4-tools/bin:/home/palak/miniconda3/bin:/home/palak/.cargo/bin:/usr/bin:/bin"
export HOME=/home/palak
# Stack down during regen -> blackout window meaningless per owner ruling 2026-10-16.
# Floors still enforced per step via courtesy.sh disk/RAM checks.
export SKIP_BLACKOUT=1
log "=== heavy regen start ==="
[ "$(id -un)" = "palak" ] || { log "FATAL: must run as palak (git ownership), got $(id -un)"; exit 8; }
[ "$(date +%F)" = "2026-10-17" ] || { log "FATAL: wrong date, aborting"; exit 8; }  # one-shot 2026-10-17 re-fire (misfire 10-16: root identity)
[ "$(df -BG --output=avail "$WT" | tail -1 | tr -dc 0-9)" -ge 150 ] || { log "FATAL: disk floor"; exit 8; }
[ "$(free -g | awk '/^Mem:/{print $7}')" -ge 12 ] || { log "FATAL: RAM floor"; exit 8; }
cd "$WT"
[ "$(git rev-parse --abbrev-ref HEAD)" = "perf/bounded-amplicon-memory" ] || { log "FATAL: wrong branch"; exit 8; }
git diff --quiet HEAD -- crates harness Cargo.toml Cargo.lock || { log "FATAL: dirty tree (fingerprint gate)"; exit 8; }
command -v conda >/dev/null && conda run -n vardict-oracle true 2>&1 | tee -a "$LOG" || { log "FATAL: oracle env"; exit 8; }
# Oracle smoke: version probes are the lightest check that the java+R+perl chain
# is alive (vardict-java launches the JVM; Rscript covers the testsomatic.R /
# teststrandbias.R stage; perl covers the var2vcf_*.pl stage). Full chain proof
# comes from the Tier-2 runs themselves below.
conda run -n vardict-oracle vardict-java -h 2>&1 | tee -a "$LOG" || { log "FATAL: oracle java probe"; exit 8; }
conda run -n vardict-oracle Rscript --version 2>&1 | tee -a "$LOG" || { log "FATAL: oracle R probe"; exit 8; }
conda run -n vardict-oracle perl --version 2>&1 | tail -n 2 | tee -a "$LOG" || { log "FATAL: oracle perl probe"; exit 8; }
rustup target list --installed 2>/dev/null | grep -q musl || { log "FATAL: musl target"; exit 8; }
# Input-asset preflight: every BAM/ref/BED/truth/sdf path referenced by the
# suite (list derived from the recorded evidence-json cmds + truth paths +
# cram_twin.sh defaults). Abort before touching the media stack.
for f in scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa \
    scratch/g4/ref/giabv3/GRCh38_GIABv3.fa \
    scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta \
    scratch/g4/amplicon/ref.fa \
    scratch/g4/ref/d1vd1/GRCh38.d1.vd1.sdf \
    scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf \
    scratch/g4/seqc2/HCC1395.bam scratch/g4/seqc2/HCC1395BL.bam \
    scratch/g4/hg008/HG008-T.slice.bam scratch/g4/hg008/HG008-N.slice.bam \
    scratch/g4/hg008/HG008-T.slice.cram scratch/g4/hg008/HG008-N.slice.cram \
    scratch/g4/hg002/HG002.chr20.bam scratch/g4/amplicon/amplicon.bam \
    scratch/g4/seqc2/exome.engine.bed scratch/g4/seqc2/slices.engine.bed \
    scratch/g4/seqc2/eval.bed scratch/g4/hg008/engine.bed scratch/g4/hg008/eval.bed \
    scratch/g4/hg002/eval.bed scratch/g4/amplicon/correct_8columns_amplicon.bed \
    scratch/g4/seqc2/truth/truth.vcf.gz \
    scratch/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz; do
  [ -e "$f" ] || { log "FATAL: missing input asset $f"; exit 8; }  # -e: SDFs are directories
done
for d in scratch/g4/seqc2/truth scratch/g4/hg008/truth; do
  [ -d "$d" ] || { log "FATAL: missing input dir $d"; exit 8; }
done
log "input-asset preflight ok"
OOMD_WAS_ACTIVE=$(systemctl is-active systemd-oomd || true)
log "oomd was: $OOMD_WAS_ACTIVE"
docker info >/dev/null || { log "FATAL: docker daemon down (down-confirm loop would be meaningless)"; exit 8; }
[ "${SKIP_STACK:-0}" = "1" ] || {
  sudo systemctl stop systemd-oomd 2>&1 | tee -a "$LOG"
  "$VENV_PY" "$LIFECYCLE" down 2>&1 | tee -a "$LOG"
  STACK_DOWN=1
  for s in deluge sonarr radarr plex tdarr sabnzbd; do
    if docker ps --format '{{.Names}}' | grep -qx "$s"; then log "FATAL: $s still Up after down"; exit 8; fi
  done
  log "media stack down confirmed"
}
assert_match() { # $1=compare.txt $2=expected_variants
  grep -q "^MATCH ($2 variants)$" "$1" || { log "FATAL: no MATCH($2) in $1"; tail -n 5 "$1" | tee -a "$LOG"; exit 9; }
}
stamp() { # $1=name $2=result(MATCH|EQUAL|PASS) $3=keepdir|$4...=record cmd
  local name=$1 res=$2 kd=$3; shift 3
  python3 harness/record_evidence.py "$name" "$res" "$kd" -- "$@" 2>&1 | tee -a "$LOG"
}
echo "--- [1] seqc2 exome Tier-2 (~70min) ---" | tee -a "$LOG"
rm -f scratch/g4/seqc2/cmp-exome/compare.txt scratch/g4/seqc2/cmp-exome/java.vcf scratch/g4/seqc2/cmp-exome/rust.vcf
EST_MIN=80 bash harness/courtesy.sh env KEEP_DIR=scratch/g4/seqc2/cmp-exome MIN_VARIANTS=2000 JAVA_OPTS=-Xmx8g \
  bash harness/compare.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" \
  scratch/g4/seqc2/exome.engine.bed "HCC1395|HCC1395BL" -th 4 2>&1 | tee -a "$LOG"
assert_match scratch/g4/seqc2/cmp-exome/compare.txt 623444
stamp seqc2-exome-tier2 MATCH scratch/g4/seqc2/cmp-exome env KEEP_DIR=scratch/g4/seqc2/cmp-exome MIN_VARIANTS=2000 JAVA_OPTS=-Xmx8g bash harness/compare.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" scratch/g4/seqc2/exome.engine.bed "HCC1395|HCC1395BL" -th 4
echo "--- [2] chr1 scale (~60min) ---" | tee -a "$LOG"
rm -f scratch/g4/scale/extra.json
EST_MIN=70 bash harness/courtesy.sh bash harness/scale_bench.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa scratch/g4/seqc2/HCC1395.bam scratch/g4/seqc2/exome.engine.bed HCC1395 scratch/g4/scale 2>&1 | tee -a "$LOG"
grep -q '"pass": true' scratch/g4/scale/extra.json || { log "FATAL: scale gate not PASS"; exit 9; }
stamp scale PASS scratch/g4/scale bash harness/scale_bench.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa scratch/g4/seqc2/HCC1395.bam scratch/g4/seqc2/exome.engine.bed HCC1395 scratch/g4/scale
echo "--- [3] hg002 Tier-2 + truth (~10min) ---" | tee -a "$LOG"
rm -f scratch/g4/hg002/cmp/compare.txt
EST_MIN=20 bash harness/courtesy.sh bash harness/run_hg002.sh 2>&1 | tee -a "$LOG"
assert_match scratch/g4/hg002/cmp/compare.txt 227562
stamp hg002-tier2 MATCH scratch/g4/hg002/cmp bash harness/compare.sh scratch/g4/ref/giabv3/GRCh38_GIABv3.fasta scratch/g4/hg002/HG002.chr20.bam scratch/g4/hg002/eval.bed HG002 -f 0.01 -th 8
echo "--- [4] hg008 Tier-2 (~10min) ---" | tee -a "$LOG"
rm -f scratch/g4/hg008/cmp/compare.txt scratch/g4/hg008/cmp/java.vcf scratch/g4/hg008/cmp/rust.vcf
EST_MIN=20 bash harness/courtesy.sh env KEEP_DIR=scratch/g4/hg008/cmp MIN_VARIANTS=1000 JAVA_OPTS=-Xmx8g \
  bash harness/compare.sh scratch/g4/ref/giabv3/GRCh38_GIABv3.fa "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam" \
  scratch/g4/hg008/engine.bed "HG008-T|HG008-N" -th 4 2>&1 | tee -a "$LOG"
assert_match scratch/g4/hg008/cmp/compare.txt 18426
stamp hg008-tier2 MATCH scratch/g4/hg008/cmp env KEEP_DIR=scratch/g4/hg008/cmp MIN_VARIANTS=1000 JAVA_OPTS=-Xmx8g bash harness/compare.sh scratch/g4/ref/giabv3/GRCh38_GIABv3.fa "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam" scratch/g4/hg008/engine.bed "HG008-T|HG008-N" -th 4
echo "--- [5] amplicon Tier-2 + truths + slices + cram + musl (~20min) ---" | tee -a "$LOG"
rm -f scratch/g4/amplicon/cmp/compare.txt
EST_MIN=10 bash harness/courtesy.sh env KEEP_DIR=scratch/g4/amplicon/cmp \
  bash harness/compare.sh scratch/g4/amplicon/ref.fa scratch/g4/amplicon/amplicon.bam scratch/g4/amplicon/correct_8columns_amplicon.bed SAMPLE -a 10:0.95 -th 8 2>&1 | tee -a "$LOG"
assert_match scratch/g4/amplicon/cmp/compare.txt 2
stamp amplicon-tier2 MATCH scratch/g4/amplicon/cmp bash harness/compare.sh scratch/g4/amplicon/ref.fa scratch/g4/amplicon/amplicon.bam scratch/g4/amplicon/correct_8columns_amplicon.bed SAMPLE -a 10:0.95 -th 8
EST_MIN=15 bash harness/courtesy.sh bash harness/concordance_pair.sh seqc2-truth scratch/g4/seqc2/cmp-exome/java.vcf scratch/g4/seqc2/cmp-exome/rust.vcf scratch/g4/seqc2/truth/truth.vcf.gz scratch/g4/seqc2/eval.bed scratch/g4/ref/d1vd1/GRCh38.d1.vd1.sdf scratch/g4/seqc2/truth HCC1395 scratch/g4/seqc2/cmp-exome/cpp.vcf 2>&1 | tee -a "$LOG"
EST_MIN=10 bash harness/courtesy.sh bash harness/concordance_pair.sh hg008-truth scratch/g4/hg008/cmp/java.vcf scratch/g4/hg008/cmp/rust.vcf scratch/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz scratch/g4/hg008/eval.bed scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf scratch/g4/hg008/truth HG008-T scratch/g4/hg008/cmp/cpp.vcf 2>&1 | tee -a "$LOG"
EST_MIN=10 bash harness/courtesy.sh bash harness/concordance_multimodal.sh 2>&1 | tee -a "$LOG"
rm -f scratch/g4/hg008/cram/compare.txt
EST_MIN=15 bash harness/courtesy.sh env KEEP_DIR=scratch/g4/hg008/cram CMP_TIMEOUT=12h bash harness/cram_twin.sh 2>&1 | tee -a "$LOG"
grep -q "^MATCH" scratch/g4/hg008/cram/compare.txt || { log "FATAL: cram twin no MATCH"; exit 9; }
stamp hg008-cram-twin MATCH scratch/g4/hg008/cram env KEEP_DIR=scratch/g4/hg008/cram CMP_TIMEOUT=12h bash harness/cram_twin.sh
cargo build --release --target x86_64-unknown-linux-musl -p vardictrust 2>&1 | tail -n 2 | tee -a "$LOG"
bash harness/bin_equiv.sh 2>&1 | tail -n 3 | tee -a "$LOG"
python3 harness/record_evidence.py musl-ci PASS - -- cargo build --release --target x86_64-unknown-linux-musl -p vardictrust 2>&1 | tee -a "$LOG"
echo "--- [6] seqc2 slices Tier-2 (fast) ---" | tee -a "$LOG"
rm -f scratch/g4/seqc2/cmp-slices/compare.txt
EST_MIN=10 bash harness/courtesy.sh env KEEP_DIR=scratch/g4/seqc2/cmp-slices MIN_VARIANTS=200 JAVA_OPTS=-Xmx8g \
  bash harness/compare.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" \
  scratch/g4/seqc2/slices.engine.bed "HCC1395|HCC1395BL" -th 8 2>&1 | tee -a "$LOG"
assert_match scratch/g4/seqc2/cmp-slices/compare.txt 4151
stamp seqc2-slices-tier2 MATCH scratch/g4/seqc2/cmp-slices env KEEP_DIR=scratch/g4/seqc2/cmp-slices MIN_VARIANTS=200 JAVA_OPTS=-Xmx8g bash harness/compare.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa "scratch/g4/seqc2/HCC1395.bam|scratch/g4/seqc2/HCC1395BL.bam" scratch/g4/seqc2/slices.engine.bed "HCC1395|HCC1395BL" -th 8
echo "--- [7] perf-budget (reads fresh tiers) ---" | tee -a "$LOG"
# perf_budget.py stdout is the JSON (must land in extra.json AND the log);
# stderr carries only diagnostics (must NOT pollute extra.json). Save stdout to
# the file, capture rc via the || guard (pipefail-safe), then mirror both to LOG.
mkdir -p scratch/g4/perf-budget
rc=0
python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 >scratch/g4/perf-budget/extra.json 2>scratch/g4/perf-budget/stderr.txt || rc=$?
cat scratch/g4/perf-budget/extra.json >>"$LOG"
cat scratch/g4/perf-budget/stderr.txt >>"$LOG"; rm -f scratch/g4/perf-budget/stderr.txt
[ "$rc" -eq 0 ] || { log "FATAL: perf_budget.py exit $rc"; exit 9; }
grep -q '"pass": true' scratch/g4/perf-budget/extra.json || { log "FATAL: perf budget not PASS"; exit 9; }
stamp perf-budget PASS scratch/g4/perf-budget python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
echo "--- [8] recertify fast suites + acceptance ---" | tee -a "$LOG"
EST_MIN=30 bash harness/courtesy.sh bash harness/recertify.sh 2>&1 | tee -a "$LOG"
log "=== heavy regen complete — review $LOG, then commit evidence + open next steps ==="
