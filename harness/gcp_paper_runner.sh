#!/usr/bin/env bash
# =============================================================================
# harness/gcp_paper_runner.sh
# Single-Instance GCP Orchestrator & Worker for `paper/gcp-study` (v1.0.3).
#
# Architecture (locked in docs/validation/gcp-study-protocol.md §6 & gcp-prereg.md §A12):
#   - Single dedicated VM: `n2-highmem-32` (32 vCPUs, 256 GB RAM, 250 GB pd-ssd boot + 8x375 GB [3.0 TB] NVMe Local SSD RAID0, us-central1-a)
#   - 2-Strike Spot -> On-Demand Escalation Policy:
#       Attempt 1: SPOT (--provisioning-model=SPOT)
#       Attempt 2: SPOT (--provisioning-model=SPOT)
#       Attempt 3: STANDARD (--provisioning-model=STANDARD, On-Demand)
#   - Per-cell GCS checkpointing (`gs://vardict-paper-trial-results/checkpoints-v103/<cell>.done`)
#     so preempted Spot runs resume cleanly without repeating completed cells.
#   - Detached systemd service on VM (`vardict-paper-worker`) + GCS heartbeat polling so
#     transient SSH drops never interrupt a healthy VM; preemption is detected via
#     `gcloud compute instances describe` status transitions.
#   - Cgroup slice isolation per cell (`systemd-run --scope -p AllowedCPUs=... -p MemoryMax=...`).
#   - `export GOLDEN_MODE=off` enforced globally on all timed/OOM cells.
# =============================================================================
set -euo pipefail

export GOLDEN_MODE=off
export CLOUDSDK_CORE_ACCOUNT="${CLOUDSDK_CORE_ACCOUNT:-vardict-agent@vardict-paper-trial.iam.gserviceaccount.com}"
export CLOUDSDK_CORE_PROJECT="${CLOUDSDK_CORE_PROJECT:-vardict-paper-trial}"

GCP_PROJECT="${CLOUDSDK_CORE_PROJECT}"
GCP_ZONE="${GCP_ZONE:-us-central1-a}"
GCP_SA="${CLOUDSDK_CORE_ACCOUNT}"
VM_NAME="${VM_NAME:-vardict-paper-n2hm32}"
MACHINE_TYPE="${MACHINE_TYPE:-n2-highmem-32}"
BOOT_DISK_SIZE="${BOOT_DISK_SIZE:-250GB}"
BOOT_DISK_TYPE="${BOOT_DISK_TYPE:-pd-ssd}"
LOCAL_SSD_COUNT="${LOCAL_SSD_COUNT:-8}"
MAX_RUN_DURATION="${MAX_RUN_DURATION:-24h}"
KEEP_VM_ON_FAIL="${KEEP_VM_ON_FAIL:-1}"

GCS_BAMS="${GCS_BAMS:-gs://vardict-paper-trial-bams}"
GCS_REFS="${GCS_REFS:-gs://vardict-paper-trial-refs}"
GCS_RESULTS="${GCS_RESULTS:-gs://vardict-paper-trial-results}"
GCS_EVIDENCE="${GCS_EVIDENCE:-gs://vardict-paper-trial-evidence}"
GCS_CKPT="${GCS_CKPT:-${GCS_RESULTS}/checkpoints-v103}"

EXPECTED_FP="${EXPECTED_FP:-27f3667b32ac14f5}"
WORK_DIR="${WORK_DIR:-/var/tmp/vardict-paper}"

usage() {
  cat <<'EOF'
Usage:
  bash harness/gcp_paper_runner.sh --plan-check
  bash harness/gcp_paper_runner.sh --controller [--max-attempts 3]
  bash harness/gcp_paper_runner.sh --worker [--stage all|1|2|3|4]

Modes:
  --plan-check   Validate 2-strike Spot->On-Demand policy, 2048GB pd-ssd config, and checkpoint logic (offline safe)
  --controller   Host-side 2-strike Spot -> On-Demand VM lifecycle manager & progress monitor
  --worker       On-VM 4-stage benchmark executor with systemd-run cgroup slices and per-cell GCS checkpoints
EOF
}

# Maps attempt number (1..3) to GCP provisioning model per Gate 8 (2-Strike Spot -> On-Demand).
provisioning_model_for_attempt() {
  local attempt="$1"
  if [ "$attempt" -le 2 ]; then
    echo "SPOT"
  else
    echo "STANDARD"
  fi
}

# Builds the exact `gcloud compute instances create` command array for a given attempt.
build_gcloud_create_cmd() {
  local attempt="$1"
  local model
  model="$(provisioning_model_for_attempt "$attempt")"
  local cmd=(
    gcloud compute instances create "$VM_NAME"
    "--project=${GCP_PROJECT}"
    "--zone=${GCP_ZONE}"
    "--machine-type=${MACHINE_TYPE}"
    "--boot-disk-size=${BOOT_DISK_SIZE}"
    "--boot-disk-type=${BOOT_DISK_TYPE}"
  )
  local i
  for ((i = 1; i <= LOCAL_SSD_COUNT; i++)); do
    cmd+=("--local-ssd=interface=NVME")
  done
  cmd+=(
    "--image-family=ubuntu-2204-lts"
    "--image-project=ubuntu-os-cloud"
    "--service-account=${GCP_SA}"
    "--scopes=cloud-platform"
    "--provisioning-model=${model}"
    "--max-run-duration=${MAX_RUN_DURATION}"
    "--instance-termination-action=DELETE"
  )
  if [ "$model" = "SPOT" ]; then
    cmd+=("--maintenance-policy=TERMINATE")
  fi
  printf '%s\n' "${cmd[@]}"
}

# Single-stream gsutil upload helper (Task 5a C1/C2: avoid gcloud storage cp hangs
# and retain MD5 on >150 MB objects via parallel_composite_upload_threshold=0).
gcs_upload() {
  local src="$1" dst="$2"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    mkdir -p "$(dirname "$dst")"
    cp -r "$src" "$dst"
  else
    gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m cp -r "$src" "$dst"
  fi
}

gcs_stat_done() {
  local cell_id="$1"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    [ -f "${LOCAL_CKPT_DIR}/${cell_id}.done" ]
  else
    gsutil -q stat "${GCS_CKPT}/${cell_id}.done" 2>/dev/null
  fi
}

gcs_mark_done() {
  local cell_id="$1"
  local stamp
  stamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    mkdir -p "${LOCAL_CKPT_DIR}"
    echo "$stamp" > "${LOCAL_CKPT_DIR}/${cell_id}.done"
  else
    printf '%s\n' "$stamp" | gsutil -q cp - "${GCS_CKPT}/${cell_id}.done"
  fi
}

gcs_write_heartbeat() {
  local stage="$1" cell_id="$2" state="$3"
  local ts
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  local payload="{\"timestamp\":\"${ts}\",\"vm\":\"${VM_NAME}\",\"stage\":\"${stage}\",\"cell\":\"${cell_id}\",\"state\":\"${state}\",\"golden_mode\":\"${GOLDEN_MODE}\"}"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    mkdir -p "${LOCAL_CKPT_DIR}"
    printf '%s\n' "$payload" > "${LOCAL_CKPT_DIR}/status.json"
  else
    printf '%s\n' "$payload" | gsutil -q cp - "${GCS_CKPT}/status.json" || true
  fi
}

# Executes a single cell inside a systemd-run CPU+RAM cgroup slice (when systemd-run is available)
# and syncs KEEP_DIR to GCS on EVERY exit path before marking <cell_id>.done on success.
run_checkpointed_cell() {
  local stage="$1" cell_id="$2" cpus="$3" mem_max="$4"
  shift 4
  if gcs_stat_done "$cell_id"; then
    echo "[SKIP] Cell '$cell_id' already checkpointed."
    return 0
  fi

  echo "=== [RUN] Stage=$stage Cell=$cell_id CPUs=0-$((cpus-1)) (${cpus}c) MemoryMax=$mem_max GOLDEN_MODE=$GOLDEN_MODE ==="
  gcs_write_heartbeat "$stage" "$cell_id" "RUNNING"

  local cell_dir="${WORK_DIR}/cells/${cell_id}"
  mkdir -p "$cell_dir"
  export KEEP_DIR="$cell_dir"
  export GOLDEN_MODE=off

  local rc=0
  if [ "${USE_SYSTEMD_SLICE:-1}" = "1" ] && command -v systemd-run >/dev/null 2>&1; then
    if [ "$(id -u)" -eq 0 ]; then
      systemd-run --scope --quiet \
        -p "AllowedCPUs=0-$((cpus - 1))" \
        -p "MemoryMax=${mem_max}" \
        -- env PATH="$PATH" HOME="$HOME" DATA_DIR="${DATA_DIR:-scratch/g4}" KEEP_DIR="$cell_dir" GOLDEN_MODE=off "$@" || rc=$?
    else
      sudo -E systemd-run --scope --quiet \
        --uid="$(id -un)" --gid="$(id -gn)" \
        -p "AllowedCPUs=0-$((cpus - 1))" \
        -p "MemoryMax=${mem_max}" \
        -- env PATH="$PATH" HOME="$HOME" DATA_DIR="${DATA_DIR:-scratch/g4}" KEEP_DIR="$cell_dir" GOLDEN_MODE=off "$@" || rc=$?
    fi
  elif command -v taskset >/dev/null 2>&1; then
    taskset -c "0-$((cpus - 1))" "$@" || rc=$?
  else
    "$@" || rc=$?
  fi

  echo "$rc" > "${cell_dir}/cell.exit_code"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    gcs_upload "$cell_dir" "${LOCAL_CKPT_DIR}/cells/${cell_id}"
  else
    gcs_upload "$cell_dir" "${GCS_RESULTS}/cells-v103/${cell_id}" || true
  fi

  if [ "$rc" -ne 0 ]; then
    gcs_write_heartbeat "$stage" "$cell_id" "FAILED_EXIT_${rc}"
    echo "[ERROR] Cell '$cell_id' exited with code $rc (logs synced)." >&2
    return "$rc"
  fi

  gcs_mark_done "$cell_id"
  gcs_write_heartbeat "$stage" "$cell_id" "DONE"
  echo "=== [DONE] Cell '$cell_id' checkpointed ==="
}

# =============================================================================
# Worker Stage Implementations (Executed on n2-highmem-32 / 2 TB pd-ssd)
# =============================================================================

ensure_submodule_excludes() {
  for ex in .git/modules/VarDictJava/info/exclude VarDictJava/.git/info/exclude; do
    if [ -f "$ex" ]; then
      grep -qF 'testdata/fastas/*.fa' "$ex" || printf '\ntestdata/fastas/*.fa\ntestdata/fastas/*.fai\n' >> "$ex"
    fi
  done
}

worker_preflight() {
  echo "=== [Worker Preflight] ==="
  [ "$GOLDEN_MODE" = "off" ] || { echo "FATAL: GOLDEN_MODE must be off (got $GOLDEN_MODE)" >&2; exit 2; }
  ensure_submodule_excludes
  mkdir -p "$WORK_DIR/data/g4" "$WORK_DIR/cells" scratch
  if [ ! -L scratch/g4 ]; then
    rm -rf scratch/g4
    ln -sfn "$WORK_DIR/data/g4" scratch/g4
  fi
  export DATA_DIR="scratch/g4"
  local avail_gb
  avail_gb="$(df -BG --output=avail "$WORK_DIR" | tail -1 | tr -dc '0-9')"
  local min_gb="${MIN_DISK_GB:-1200}"
  if [ "$avail_gb" -lt "$min_gb" ]; then
    echo "FATAL: $WORK_DIR has ${avail_gb}GB free, need >= ${min_gb}GB (NVMe RAID0 SSD required)" >&2
    exit 2
  fi
  local fp
  fp="$(bash harness/engine_fingerprint.sh)"
  if [ "$fp" != "$EXPECTED_FP" ]; then
    echo "FATAL: engine fingerprint mismatch: got '$fp', expected '$EXPECTED_FP'" >&2
    exit 2
  fi
  cargo build --release -p vardictrust
}

run_stage_1() {
  echo "=== [Stage 1: Task 5b Pre-Full-Genome Gates] ==="
  export DATA_DIR="scratch/g4"

  run_checkpointed_cell "stage1" "s1_stage_assets" 32 "128G" bash -c "
    set -euo pipefail
    if gsutil -q stat '${GCS_BAMS}/g4-cache/s1_ready.done' 2>/dev/null; then
      gsutil -m rsync -r '${GCS_BAMS}/g4-cache/amplicon' scratch/g4/amplicon
      gsutil -m rsync -r '${GCS_BAMS}/g4-cache/hg002' scratch/g4/hg002
      gsutil -m rsync -r '${GCS_BAMS}/g4-cache/hg008' scratch/g4/hg008
    fi
    bash harness/stage_g4.sh refs amplicon hg002 hg008
    gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m rsync -r scratch/g4/amplicon '${GCS_BAMS}/g4-cache/amplicon' || true
    gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m rsync -r scratch/g4/hg002 '${GCS_BAMS}/g4-cache/hg002' || true
    gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m rsync -r scratch/g4/hg008 '${GCS_BAMS}/g4-cache/hg008' || true
    date -u +%Y-%m-%dT%H:%M:%SZ | gsutil -q cp - '${GCS_BAMS}/g4-cache/s1_ready.done' || true
  "
  ensure_submodule_excludes

  # H1a Simple-mode tile-size equivalence (500k, 1M, 2M, 10M vs untiled oracle on HG002 slice sub-window)
  run_checkpointed_cell "stage1" "s1_h1a_simple_tiles" 8 "32G" bash -c '
    set -euo pipefail
    head -n 40 "$DATA_DIR/hg002/eval.bed" > "$KEEP_DIR/h1a_simple.bed"
    ORACLE_XMX_GB=22 bash harness/run_oracle.sh \
      "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" \
      "$DATA_DIR/hg002/HG002.chr20.bam" \
      "$KEEP_DIR/h1a_simple.bed" \
      "$KEEP_DIR/oracle_untiled.vcf" HG002 -f 0.01 -th 8 -z
    for ts in 500000 1000000 2000000 10000000; do
      target/release/vardictrust \
        -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" \
        -f 0.01 -N HG002 \
        -b "$DATA_DIR/hg002/HG002.chr20.bam" \
        -c 1 -S 2 -E 3 -g 4 -z -th 8 \
        --tile-size "$ts" \
        "$KEEP_DIR/h1a_simple.bed" > "$KEEP_DIR/rust_tile_${ts}.vcf"
      python3 harness/compare_vcf.py --min-variants 100 \
        "$KEEP_DIR/oracle_untiled.vcf" "$KEEP_DIR/rust_tile_${ts}.vcf" \
        | tee "$KEEP_DIR/compare_tile_${ts}.txt"
      grep -q "^MATCH" "$KEEP_DIR/compare_tile_${ts}.txt"
    done
  '

  # Gate 7: Multi-region reference offset check (2+ distinct regions in one invocation)
  run_checkpointed_cell "stage1" "s1_gate7_ref_multiregion" 8 "32G" bash -c '
    set -euo pipefail
    sed -n "1,15p;50,65p" "$DATA_DIR/hg002/eval.bed" > "$KEEP_DIR/gate7_multiregion.bed"
    GOLDEN_MODE=off MIN_VARIANTS=50 ORACLE_XMX_GB=22 \
      bash harness/compare.sh \
      "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" \
      "$DATA_DIR/hg002/HG002.chr20.bam" \
      "$KEEP_DIR/gate7_multiregion.bed" HG002 -f 0.01 -th 8
  '

  # Gate 1 Ceiling Cell A: Amplicon deep panel (-th 8, assert Rust peak RSS <= 8 GB)
  run_checkpointed_cell "stage1" "s1_ceiling_amplicon" 8 "32G" bash -c '
    set -euo pipefail
    GOLDEN_MODE=off ORACLE_XMX_GB=22 \
      bash harness/compare.sh \
      "$DATA_DIR/amplicon/ref.fa" \
      "$DATA_DIR/amplicon/amplicon.bam" \
      "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" SAMPLE -a 10:0.95 -th 8
    rss_kb=$(awk -F": " "/Maximum resident set size/{print \$2}" "$KEEP_DIR/rust.time" | tr -dc "0-9")
    echo "amplicon_ceiling_rss_kb=$rss_kb" | tee "$KEEP_DIR/ceiling_amplicon.txt"
    [ "$rss_kb" -le 8000000 ] || { echo "FATAL: amplicon ceiling RSS ${rss_kb} KB > 8 GB" >&2; exit 9; }
  '

  # Gate 1 Ceiling Cell B: HG008 195x/169x somatic WGS slice (-th 8, assert Rust peak RSS <= 8 GB)
  run_checkpointed_cell "stage1" "s1_ceiling_hg008" 8 "32G" bash -c '
    set -euo pipefail
    GOLDEN_MODE=off MIN_VARIANTS=1000 ORACLE_XMX_GB=22 \
      bash harness/compare.sh \
      "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" \
      "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" \
      "$DATA_DIR/hg008/engine.bed" "HG008-T|HG008-N" -th 8
    rss_kb=$(awk -F": " "/Maximum resident set size/{print \$2}" "$KEEP_DIR/rust.time" | tr -dc "0-9")
    echo "hg008_ceiling_rss_kb=$rss_kb" | tee "$KEEP_DIR/ceiling_hg008.txt"
    [ "$rss_kb" -le 8000000 ] || { echo "FATAL: HG008 ceiling RSS ${rss_kb} KB > 8 GB" >&2; exit 9; }
  '
}

run_stage_2() {
  echo "=== [Stage 2: Task 5c Tier-P 5-Sample Replication + Fig-2 Truth + Table-3 -R chr1] ==="
  local tp_dir="${WORK_DIR}/data/tier-p"
  mkdir -p "$tp_dir/bams" "$tp_dir/truth" "$tp_dir/archive"

  run_checkpointed_cell "stage2" "s2_stage_tierp_assets" 32 "128G" bash -c "
    set -euo pipefail
    gsutil -m cp '${GCS_BAMS}/tier-p/*' '$tp_dir/bams/'
    gsutil -m cp '${GCS_REFS}/tier-p-truth/*' '$tp_dir/truth/'
    if [ ! -f '$tp_dir/archive/Archive.zip' ]; then
      gsutil cp '${GCS_REFS}/vardictcpp/Archive.zip' '$tp_dir/archive/Archive.zip'
    fi
    unzip -q -o '$tp_dir/archive/Archive.zip' 'reference/*' 'beds/*' 'scripts/*' 'performance_benchmark/*' 'MD5SUMS.txt' -d '$tp_dir/archive/'
    if [ ! -x '$tp_dir/vardictcpp-ae93fe39dc/build/src/vardictcpp' ]; then
      rm -rf '$tp_dir/vardictcpp-ae93fe39dc'
      git clone https://github.com/MHH-Bioinformatics-Hematology/vardictcpp.git '$tp_dir/vardictcpp-ae93fe39dc'
      cd '$tp_dir/vardictcpp-ae93fe39dc'
      git checkout ae93fe39dc
      conda run -n vardict-oracle bash -c 'cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 16'
    fi
  "

  # Derive covered BEDs via mkbed() and run 3-way Table-1 parity + Table-2 speed across 5 SRR samples
  for acc in SRR15006386 SRR15006375 SRR15006376 SRR15006540 SRR8657348; do
    run_checkpointed_cell "stage2" "s2_tierp_3way_${acc}" 8 "32G" bash -c "
      set -euo pipefail
      REF='$tp_dir/archive/reference/hg19.fa'
      BAM='$tp_dir/bams/${acc}.bam'
      BED='$tp_dir/bams/${acc}.covered.bed'
      if [ ! -s \"\$BED\" ]; then
        bedtools bamtobed -i \"\$BAM\" | awk '\$5 >= 1' | sort -k1,1 -k2,2n | bedtools merge -d 100 -i - | awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\"R\"}' > \"\$BED\"
      fi
      wc -l \"\$BED\" > \"\$KEEP_DIR/${acc}.bed_lines\"
      for th in 1 8; do
        for run in 1 2 3; do
          /usr/bin/time -v -o \"\$KEEP_DIR/rust_th\${th}_run\${run}.time\" \
            target/release/vardictrust -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -th \"\$th\" \"\$BED\" > \"\$KEEP_DIR/rust_th\${th}.vcf\"
          /usr/bin/time -v -o \"\$KEEP_DIR/cpp_th\${th}_run\${run}.time\" \
            '$tp_dir/vardictcpp-ae93fe39dc/build/src/vardictcpp' -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z --chunk 50000 -th \"\$th\" \"\$BED\" -vcf > \"\$KEEP_DIR/cpp_th\${th}.vcf\"
          ORACLE_XMX_GB=22 /usr/bin/time -v -o \"\$KEEP_DIR/java_th\${th}_run\${run}.time\" \
            bash harness/run_oracle.sh \"\$REF\" \"\$BAM\" \"\$BED\" \"\$KEEP_DIR/java_th\${th}.vcf\" '${acc}' -f 0.01 -z -th \"\$th\"
        done
        python3 harness/compare_vcf.py \"\$KEEP_DIR/java_th\${th}.vcf\" \"\$KEEP_DIR/rust_th\${th}.vcf\" | tee \"\$KEEP_DIR/compare_th\${th}.txt\"
        grep -q '^MATCH' \"\$KEEP_DIR/compare_th\${th}.txt\"
      done
    "
  done

  # Somatic pair SRR15006376 (tumor) | SRR15006386 (normal)
  run_checkpointed_cell "stage2" "s2_tierp_somatic_76_86" 8 "32G" bash -c "
    set -euo pipefail
    REF='$tp_dir/archive/reference/hg19.fa'
    TUM='$tp_dir/bams/SRR15006376.bam'
    NOR='$tp_dir/bams/SRR15006386.bam'
    BED='$tp_dir/bams/SRR15006376.covered.bed'
    for th in 1 8; do
      for run in 1 2 3; do
        /usr/bin/time -v -o \"\$KEEP_DIR/rust_som_th\${th}_run\${run}.time\" \
          target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 -z -th \"\$th\" \"\$BED\" > \"\$KEEP_DIR/rust_som_th\${th}.vcf\"
        /usr/bin/time -v -o \"\$KEEP_DIR/cpp_som_th\${th}_run\${run}.time\" \
          '$tp_dir/vardictcpp-ae93fe39dc/build/src/vardictcpp' -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 -z --chunk 50000 -th \"\$th\" \"\$BED\" -vcf > \"\$KEEP_DIR/cpp_som_th\${th}.vcf\"
        ORACLE_XMX_GB=22 /usr/bin/time -v -o \"\$KEEP_DIR/java_som_th\${th}_run\${run}.time\" \
          bash harness/run_oracle.sh \"\$REF\" \"\${TUM}|\${NOR}\" \"\$BED\" \"\$KEEP_DIR/java_som_th\${th}.vcf\" 'SRR15006376|SRR15006386' -f 0.01 -z -th \"\$th\"
      done
      python3 harness/compare_vcf.py \"\$KEEP_DIR/java_som_th\${th}.vcf\" \"\$KEEP_DIR/rust_som_th\${th}.vcf\" | tee \"\$KEEP_DIR/compare_som_th\${th}.txt\"
      grep -q '^MATCH' \"\$KEEP_DIR/compare_som_th\${th}.txt\"
    done
  "

  # Wolff 2026 Table-3 Direct -R chr1:1-249250621 Whole-Chromosome Blowup Cell on SRR8657348
  run_checkpointed_cell "stage2" "s2_table3_R_chr1_rust_and_cpp" 8 "32G" bash -c "
    set -euo pipefail
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    for run in 1 2 3; do
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_R_chr1_run\${run}.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 > \"\$KEEP_DIR/rust_R_chr1.vcf\"
      rc=0
      /usr/bin/time -v -o \"\$KEEP_DIR/cpp_R_chr1_run\${run}.time\" \
        '$tp_dir/vardictcpp-ae93fe39dc/build/src/vardictcpp' -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 -vcf > \"\$KEEP_DIR/cpp_R_chr1.vcf\" || rc=\$?
      echo \"\$rc\" > \"\$KEEP_DIR/cpp_R_chr1_run\${run}.exit\"
    done
  "

  # Table-3 Java -Xmx4g OOM reproduction & -Xmx200g completion cell
  run_checkpointed_cell "stage2" "s2_table3_R_chr1_java_200g" 8 "240G" bash -c "
    set -euo pipefail
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    rc4=0
    _JAVA_OPTIONS='-Xmx4g' /usr/bin/time -v -o \"\$KEEP_DIR/java_xmx4g_R_chr1.time\" \
      conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 > /dev/null 2>\"\$KEEP_DIR/java_xmx4g.err\" || rc4=\$?
    echo \"\$rc4\" > \"\$KEEP_DIR/java_xmx4g.exit\"
    _JAVA_OPTIONS='-Xmx200g' /usr/bin/time -v -o \"\$KEEP_DIR/java_xmx200g_R_chr1.time\" \
      conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 > \"\$KEEP_DIR/java_xmx200g_R_chr1.var\" 2>\"\$KEEP_DIR/java_xmx200g.err\"
  "
}

run_stage_3() {
  echo "=== [Stage 3: Task 5d Cross-Tier Full-Sample RNA Arm (HG002 GM24385 51.8 GiB BAM)] ==="
  export DATA_DIR="scratch/g4"
  local rna_dir="${WORK_DIR}/data/rna"
  mkdir -p "$rna_dir"

  run_checkpointed_cell "stage3" "s3_stage_rna_aligned_bam" 32 "128G" bash -c "
    set -euo pipefail
    if [ ! -f '$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam' ]; then
      gsutil -m cp '${GCS_BAMS}/rna-hg002/aligned/HG002_rep1_Aligned.sortedByCoord.out.bam*' '$rna_dir/'
    fi
    samtools quickcheck '$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam'
  "

  run_checkpointed_cell "stage3" "s3_rna_splice_parity" 8 "64G" bash -c "
    set -euo pipefail
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    BAM='$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam'
    BED=\"\$DATA_DIR/hg002/eval.bed\"
    /usr/bin/time -v -o \"\$KEEP_DIR/rust_splice.time\" \
      target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002_RNA -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -i --emit-intermediate -th 8 \"\$BED\" > \"\$KEEP_DIR/rust.splice.tsv\"
    _JAVA_OPTIONS='-Xmx44g' /usr/bin/time -v -o \"\$KEEP_DIR/java_splice.time\" \
      conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N HG002_RNA -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z 1 -i -th 8 \"\$BED\" > \"\$KEEP_DIR/java.splice.tsv\"
    LC_ALL=C sort \"\$KEEP_DIR/java.splice.tsv\" > \"\$KEEP_DIR/java.splice.sorted.tsv\"
    LC_ALL=C sort \"\$KEEP_DIR/rust.splice.tsv\" > \"\$KEEP_DIR/rust.splice.sorted.tsv\"
    diff -u \"\$KEEP_DIR/java.splice.sorted.tsv\" \"\$KEEP_DIR/rust.splice.sorted.tsv\" > \"\$KEEP_DIR/splice_diff.txt\"
    wc -l \"\$KEEP_DIR/rust.splice.sorted.tsv\" > \"\$KEEP_DIR/splice_count.txt\"
  "
}

run_stage_4() {
  echo "=== [Stage 4: Task 5e Tier-E Full SEQC2 Exome & Tier-G Full WGS (1..32 Threads)] ==="
  export DATA_DIR="scratch/g4"

  run_checkpointed_cell "stage4" "s4_stage_seqc2_full" 32 "128G" bash -c "
    set -euo pipefail
    if gsutil -q stat '${GCS_BAMS}/g4-cache/s4_ready.done' 2>/dev/null; then
      mkdir -p scratch/g4/seqc2 scratch/g4/multimodal
      gsutil -m rsync -r '${GCS_BAMS}/g4-cache/seqc2' scratch/g4/seqc2
      gsutil -m rsync -r '${GCS_BAMS}/g4-cache/multimodal' scratch/g4/multimodal
    fi
    bash harness/stage_g4.sh seqc2 multimodal
    gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m rsync -r -x 'fastq/.*' scratch/g4/seqc2 '${GCS_BAMS}/g4-cache/seqc2' || true
    gsutil -o GSUtil:parallel_composite_upload_threshold=0 -m rsync -r scratch/g4/multimodal '${GCS_BAMS}/g4-cache/multimodal' || true
    date -u +%Y-%m-%dT%H:%M:%SZ | gsutil -q cp - '${GCS_BAMS}/g4-cache/s4_ready.done' || true
  "
  ensure_submodule_excludes

  # Tier-E SEQC2 full somatic exome 3-way parity + truth
  run_checkpointed_cell "stage4" "s4_tier_e_seqc2_exome_3way" 8 "32G" bash -c '
    set -euo pipefail
    GOLDEN_MODE=off MIN_VARIANTS=2000 ORACLE_XMX_GB=22 \
      bash harness/compare.sh \
      "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa" \
      "$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam" \
      "$DATA_DIR/seqc2/exome.engine.bed" "HCC1395|HCC1395BL" -th 8
    grep -q "^MATCH (623444 variants)$" "$KEEP_DIR/compare.txt"
  '

  # Multi-core scaling 1, 2, 4, 8, 16, 32 threads on identical n2-highmem-32 silicon
  for th in 1 2 4 8 16 32; do
    local mem_cap="32G"
    [ "$th" -ge 16 ] && mem_cap="64G"
    [ "$th" -ge 32 ] && mem_cap="128G"
    run_checkpointed_cell "stage4" "s4_scale_exome_th${th}" "$th" "$mem_cap" bash -c "
      set -euo pipefail
      for run in 1 2 3; do
        /usr/bin/time -v -o \"\$KEEP_DIR/rust_th${th}_run\${run}.time\" \
          target/release/vardictrust \
          -G \"\$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa\" \
          -f 0.01 -N HCC1395 \
          -b \"\$DATA_DIR/seqc2/HCC1395.bam\" \
          -c 1 -S 2 -E 3 -g 4 -z -th ${th} \
          \"\$DATA_DIR/seqc2/exome.engine.bed\" > \"\$KEEP_DIR/rust_th${th}.vcf\"
      done
    "
  done

  # Tier-G chr1 span scaling (1 Mbp, 10 Mbp, 100 Mbp, 248,956,422 bp) at 8, 16, 32 threads
  for th in 8 16 32; do
    local mem_cap="64G"
    [ "$th" -ge 32 ] && mem_cap="128G"
    run_checkpointed_cell "stage4" "s4_span_chr1_th${th}" "$th" "$mem_cap" bash -c "
      set -euo pipefail
      REF=\"\$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa\"
      BAM=\"\$DATA_DIR/seqc2/HCC1395.bam\"
      for sz in 1000000 10000000 100000000 248956422; do
        for run in 1 2 3; do
          /usr/bin/time -v -o \"\$KEEP_DIR/rust_R\${sz}_th${th}_run\${run}.time\" \
            target/release/vardictrust -G \"\$REF\" -f 0.01 -N HCC1395 -b \"\$BAM\" -th ${th} -R \"chr1:1-\${sz}\" > \"\$KEEP_DIR/rust_R\${sz}_th${th}.vcf\"
        done
      done
    "
  done

  # Concurrent jobs/node cell: 2 x 4-thread jobs inside an n2-standard-8-equivalent (8c / 32G) slice
  run_checkpointed_cell "stage4" "s4_jobs_per_node_2x4c" 8 "32G" bash -c '
    set -euo pipefail
    REF="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    BAM="$DATA_DIR/seqc2/HCC1395.bam"
    BED="$DATA_DIR/seqc2/exome.engine.bed"
    taskset -c 0-3 /usr/bin/time -v -o "$KEEP_DIR/job1_4c.time" \
      target/release/vardictrust -G "$REF" -f 0.01 -N HCC1395_J1 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th 4 "$BED" > "$KEEP_DIR/job1.vcf" &
    p1=$!
    taskset -c 4-7 /usr/bin/time -v -o "$KEEP_DIR/job2_4c.time" \
      target/release/vardictrust -G "$REF" -f 0.01 -N HCC1395_J2 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th 4 "$BED" > "$KEEP_DIR/job2.vcf" &
    p2=$!
    wait "$p1" "$p2"
  '

  # Full multimodal + extended coverage + CRAM twin recertification on the GCP host
  run_checkpointed_cell "stage4" "s4_multimodal_and_cram_twin" 8 "32G" bash -c '
    set -euo pipefail
    bash harness/run_extended_coverage.sh
    bash harness/run_multimodal_slices.sh
    CMP_TIMEOUT=12h bash harness/cram_twin.sh
    bash harness/concordance_multimodal.sh
  '
}

run_worker() {
  local target_stage="${1:-all}"
  worker_preflight
  case "$target_stage" in
    1) run_stage_1 ;;
    2) run_stage_2 ;;
    3) run_stage_3 ;;
    4) run_stage_4 ;;
    all)
      run_stage_1
      run_stage_2
      run_stage_3
      run_stage_4
      gcs_mark_done "ALL_STAGES"
      gcs_write_heartbeat "complete" "ALL_STAGES" "COMPLETED"
      ;;
    *)
      echo "Unknown worker stage: $target_stage" >&2
      exit 2
      ;;
  esac
}

# =============================================================================
# Host-Side Controller (2-Strike Spot -> On-Demand Escalation)
# =============================================================================

cleanup_vm() {
  if gcloud compute instances describe "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" >/dev/null 2>&1; then
    echo "[Controller] Deleting VM '$VM_NAME' in '$GCP_ZONE'..."
    gcloud compute instances delete "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" --quiet || true
  fi
}

vm_status() {
  gcloud compute instances describe "$VM_NAME" \
    --project="$GCP_PROJECT" \
    --zone="$GCP_ZONE" \
    --format="value(status)" 2>/dev/null || echo "ABSENT"
}

vm_provisioning_model() {
  gcloud compute instances describe "$VM_NAME" \
    --project="$GCP_PROJECT" \
    --zone="$GCP_ZONE" \
    --format="value(scheduling.provisioningModel)" 2>/dev/null || echo "ABSENT"
}

run_controller() {
  local max_attempts="${1:-3}"

  for ((attempt = 1; attempt <= max_attempts; attempt++)); do
    if gcs_stat_done "ALL_STAGES"; then
      echo "[Controller] ALL_STAGES.done already present in ${GCS_CKPT}. Deleting VM and finishing."
      cleanup_vm
      return 0
    fi

    local model
    model="$(provisioning_model_for_attempt "$attempt")"
    echo "================================================================================"
    echo " [Controller] Attempt ${attempt}/${max_attempts}: ${VM_NAME} (${MACHINE_TYPE}, ${BOOT_DISK_SIZE} ${BOOT_DISK_TYPE} + ${LOCAL_SSD_COUNT}x375GB NVMe, model=${model})"
    echo "================================================================================"

    local cur_status cur_model
    cur_status="$(vm_status)"
    cur_model="$(vm_provisioning_model)"
    if [ "$cur_status" = "RUNNING" ] && [ "$cur_model" = "$model" ]; then
      echo "[Controller] Reusing already-running VM '${VM_NAME}' (status=${cur_status}, model=${cur_model})"
    else
      cleanup_vm
      mapfile -t create_cmd < <(build_gcloud_create_cmd "$attempt")
      "${create_cmd[@]}"
    fi

    # Wait for SSH readiness
    local ssh_ok=0
    for _ in $(seq 1 36); do
      if gcloud compute ssh "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" --command="true" --quiet 2>/dev/null; then
        ssh_ok=1
        break
      fi
      sleep 10
    done
    if [ "$ssh_ok" -ne 1 ]; then
      echo "[Controller] SSH failed to come up on attempt $attempt" >&2
      continue
    fi

    # Format/mount 8x375GB NVMe Local SSD RAID0 at $WORK_DIR and start detached systemd worker
    echo "[Controller] Mounting 3.0 TB NVMe RAID0 at $WORK_DIR and starting vardict-paper-worker.service..."
    gcloud compute ssh "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" --quiet --command="
      set -euo pipefail
      sudo mkdir -p '$WORK_DIR'
      if ! mountpoint -q '$WORK_DIR'; then
        if ! command -v mdadm >/dev/null 2>&1; then
          sudo apt-get update -qq
          sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq mdadm
        fi
        if [ ! -e /dev/md0 ]; then
          sudo mdadm --create /dev/md0 --level=0 --raid-devices=${LOCAL_SSD_COUNT} /dev/disk/by-id/google-local-nvme-ssd-* --run --force
          sudo mkfs.ext4 -F -m 0 /dev/md0
        fi
        sudo mount -o discard,defaults /dev/md0 '$WORK_DIR'
      fi
      sudo chown -R \"\$(id -un):\$(id -gn)\" '$WORK_DIR'
      cat > '$WORK_DIR/bootstrap_and_run.sh' <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
export WORK_DIR='${WORK_DIR}'
export GCS_CKPT='${GCS_CKPT}'
export GOLDEN_MODE=off
exec > >(tee -a \"\${WORK_DIR}/worker.log\") 2>&1
echo \"[VM Bootstrap] Starting at \$(date -u +%Y-%m-%dT%H:%M:%SZ)\"
if ! command -v git >/dev/null || ! command -v gsutil >/dev/null || ! command -v cmake >/dev/null; then
  sudo apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq git curl unzip build-essential cmake pkg-config zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev time util-linux
fi
if [ ! -d \"\$HOME/miniconda3\" ]; then
  curl -sSLf https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
  bash /tmp/miniconda.sh -b -p \"\$HOME/miniconda3\"
fi
export PATH=\"\$HOME/miniconda3/bin:\$HOME/.cargo/bin:\$PATH\"
conda config --set auto_accept_tos yes >/dev/null 2>&1 || true
if ! command -v cargo >/dev/null; then
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.96.0
  export PATH=\"\$HOME/.cargo/bin:\$PATH\"
fi
if [ ! -d \"\${WORK_DIR}/VarDictRust\" ]; then
  git clone --branch paper/gcp-study --recurse-submodules https://github.com/palakpsheth/VarDictRust.git \"\${WORK_DIR}/VarDictRust\"
fi
cd \"\${WORK_DIR}/VarDictRust\"
git fetch origin paper/gcp-study && git checkout paper/gcp-study && git reset --hard origin/paper/gcp-study
git submodule update --init --recursive
if ! conda env list | grep -q '^vardict-oracle '; then
  conda env create -f benchmarks/env/vardict-oracle.yml -y
fi
if ! conda env list | grep -q '^g4-tools '; then
  conda env create -f benchmarks/env/g4-tools.yml -y
fi
export PATH=\"\$HOME/miniconda3/envs/g4-tools/bin:\$HOME/miniconda3/bin:\$HOME/.cargo/bin:\$PATH\"
bash harness/gcp_paper_runner.sh --worker all
EOS
      chmod +x '$WORK_DIR/bootstrap_and_run.sh'
      sudo systemctl stop vardict-paper-worker.service 2>/dev/null || true
      sudo systemctl reset-failed vardict-paper-worker.service 2>/dev/null || true
      sudo systemd-run --unit=vardict-paper-worker \
        --uid=\"\$(id -un)\" --gid=\"\$(id -gn)\" \
        --setenv=HOME=\"\$HOME\" \
        bash '$WORK_DIR/bootstrap_and_run.sh'
    "

    # Poll VM status & GCS heartbeat every 30s until ALL_STAGES.done, worker failure, or Spot preemption
    while true; do
      sleep 30
      if gcs_stat_done "ALL_STAGES"; then
        echo "[Controller] Worker completed ALL_STAGES on attempt ${attempt} (${model})!"
        cleanup_vm
        return 0
      fi
      local st
      st="$(vm_status)"
      if [ "$st" != "RUNNING" ]; then
        echo "[Controller] VM status is '${st}' on attempt ${attempt} (${model}) — Spot preemption or termination detected."
        break
      fi
      # Check if systemd unit is still active on the VM
      local unit_state
      unit_state="$(gcloud compute ssh "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" --quiet --command="systemctl is-active vardict-paper-worker.service || true" 2>/dev/null || echo "ssh_transient")"
      local hb
      hb="$(gsutil cat "${GCS_CKPT}/status.json" 2>/dev/null || echo "{}")"
      echo "[Controller $(date -u +%H:%M:%SZ)] attempt=${attempt}(${model}) vm=${st} unit=${unit_state} heartbeat=${hb}"
      if [ "$unit_state" = "failed" ] || [ "$unit_state" = "inactive" ]; then
        if gcs_stat_done "ALL_STAGES"; then
          echo "[Controller] ALL_STAGES completed!"
          cleanup_vm
          return 0
        fi
        echo "[Controller] Worker unit exited (${unit_state}) without ALL_STAGES.done; fetching tail of worker.log:" >&2
        gcloud compute ssh "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" --quiet --command="tail -n 100 '$WORK_DIR/worker.log' || true" >&2 || true
        if [ "$KEEP_VM_ON_FAIL" != "1" ]; then
          cleanup_vm
        else
          echo "[Controller] Keeping VM '${VM_NAME}' alive for inspection (KEEP_VM_ON_FAIL=1)." >&2
        fi
        return 1
      fi
    done
  done

  echo "[Controller] Exhausted ${max_attempts} attempts without completing ALL_STAGES." >&2
  cleanup_vm
  return 1
}

run_plan_check() {
  echo "=== [gcp_paper_runner: --plan-check] ==="
  [ "$GOLDEN_MODE" = "off" ] || { echo "FAIL: GOLDEN_MODE != off" >&2; exit 1; }
  [ "$(provisioning_model_for_attempt 1)" = "SPOT" ] || { echo "FAIL: attempt 1 must be SPOT" >&2; exit 1; }
  [ "$(provisioning_model_for_attempt 2)" = "SPOT" ] || { echo "FAIL: attempt 2 must be SPOT" >&2; exit 1; }
  [ "$(provisioning_model_for_attempt 3)" = "STANDARD" ] || { echo "FAIL: attempt 3 must escalate to STANDARD" >&2; exit 1; }

  local cmd1 cmd3
  cmd1="$(build_gcloud_create_cmd 1 | tr '\n' ' ')"
  cmd3="$(build_gcloud_create_cmd 3 | tr '\n' ' ')"
  [[ "$cmd1" == *"--machine-type=n2-highmem-32"* ]] || { echo "FAIL: machine-type not n2-highmem-32" >&2; exit 1; }
  [[ "$cmd1" == *"--boot-disk-size=250GB"* ]] || { echo "FAIL: boot-disk-size not 250GB" >&2; exit 1; }
  [[ "$cmd1" == *"--boot-disk-type=pd-ssd"* ]] || { echo "FAIL: boot-disk-type not pd-ssd" >&2; exit 1; }
  [[ "$cmd1" == *"--local-ssd=interface=NVME"* ]] || { echo "FAIL: missing --local-ssd=interface=NVME" >&2; exit 1; }
  [[ "$cmd1" == *"--provisioning-model=SPOT"* ]] || { echo "FAIL: attempt 1 cmd missing SPOT" >&2; exit 1; }
  [[ "$cmd3" == *"--provisioning-model=STANDARD"* ]] || { echo "FAIL: attempt 3 cmd missing STANDARD" >&2; exit 1; }
  [[ "$cmd3" != *"--maintenance-policy=TERMINATE"* ]] || { echo "FAIL: STANDARD cmd must not set maintenance-policy=TERMINATE" >&2; exit 1; }

  local tmp
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' RETURN
  export LOCAL_CKPT_DIR="$tmp/ckpt"
  export WORK_DIR="$tmp/work"
  export USE_SYSTEMD_SLICE=0

  # Verify cell execution + checkpoint creation + idempotent skip on second call
  local count_file="$tmp/run_count"
  echo 0 > "$count_file"
  run_checkpointed_cell "test_stage" "cell_alpha" 2 "4G" bash -c "echo \$(( \$(cat '$count_file') + 1 )) > '$count_file'"
  [ "$(cat "$count_file")" -eq 1 ] || { echo "FAIL: cell_alpha did not run" >&2; exit 1; }
  [ -f "$LOCAL_CKPT_DIR/cell_alpha.done" ] || { echo "FAIL: cell_alpha.done missing" >&2; exit 1; }

  run_checkpointed_cell "test_stage" "cell_alpha" 2 "4G" bash -c "echo \$(( \$(cat '$count_file') + 1 )) > '$count_file'"
  [ "$(cat "$count_file")" -eq 1 ] || { echo "FAIL: checkpointed cell_alpha ran a second time" >&2; exit 1; }

  echo "PASS gcp_paper_runner --plan-check"
}

case "${1:-}" in
  --plan-check)
    run_plan_check
    ;;
  --controller)
    shift
    max_att=3
    if [ "${1:-}" = "--max-attempts" ]; then max_att="$2"; shift 2; fi
    run_controller "$max_att"
    ;;
  --worker)
    shift
    stage="all"
    if [ "${1:-}" = "--stage" ]; then stage="$2"; shift 2; elif [ -n "${1:-}" ]; then stage="$1"; shift; fi
    run_worker "$stage"
    ;;
  -h|--help|"")
    usage
    [ -n "${1:-}" ] && exit 0 || exit 1
    ;;
  *)
    echo "Unknown option: $1" >&2
    usage >&2
    exit 1
    ;;
esac
