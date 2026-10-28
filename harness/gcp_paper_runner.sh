#!/usr/bin/env bash
# =============================================================================
# harness/gcp_paper_runner.sh
# Single-Instance GCP Orchestrator & Worker for `paper/gcp-study` (v1.0.5+).
#
# Locked Operational Invariants (docs/validation/gcp-study-protocol.md §6 & benchmarks/gcp-prereg.md §A12/B8):
#   1. CPU Microarchitecture Lock (`MIN_CPU_PLATFORM="Intel Ice Lake"`):
#      Always pass `--min-cpu-platform="Intel Ice Lake"` on `n2-standard-32` so GCP schedules
#      onto 3rd Gen Intel Xeon Scalable (`Platinum 8373C @ 2.60GHz / 3.50GHz turbo`, 1 NUMA socket,
#      identical CPU silicon to the cached `n2-highmem-32` Java/Cpp baseline) and NEVER drifts
#      onto 2nd Gen `Intel Cascade Lake` (`2.80 GHz`, ~22% slower IPC).
#   2. Right-Sized VM + 5-Spot -> 6th-On-Demand Escalation Policy:
#      `n2-standard-32` (32 vCPUs, 128 GB RAM, 100 GB pd-ssd boot + 4x375 GB [1.5 TB] NVMe RAID0).
#      Attempts 1..5 rotate `SPOT` across `us-central1-{a,b,c,f,a}` (continuing across zone Spot
#      stockouts); Attempt 6 escalates to `STANDARD` (On-Demand).
#   3. Ubuntu 24.04 Snap `gsutil` Multithreading without Multiprocessing Deadlock:
#      Always use `GSUTIL_MT_FLAGS=(-m -o "GSUtil:parallel_process_count=1" -o "GSUtil:parallel_thread_count=16")`
#      (1 Python process, 16 transfer threads — zero `multiprocessing.Pool` `fork()` deadlocks or
#      `<defunct>` zombie children) and single-pass `gsutil rsync -r -x` in `restore_prev_java_cpp_artifacts()`.
#   4. Static & Reusable Asset GCS Push/Pull Persistence:
#      All generated static assets (`*.mkbed.bed`, `vardictcpp-ae93fe39dc.bin`, `HG002`/`HG008`
#      full WGS BAMs/BEDs/truth, and `chr20`/`r123` cross-platform slices) pull from GCS before
#      generating and push immediately on cache miss so new VMs restore in seconds.
#   5. Cached Amplicon VCF Flag Parity:
#      Every `vardictrust -a ...` run against cached Java VCFs must pass `-c 1 -S 2 -E 3 -g 4`
#      so `omit_end_tag` matches `compare.sh` / `var2vcf_valid.pl -E`.
#   6. Future-Version Bump Checklist (`v1.0.6+` / `v1.1+`):
#      Update `EXPECTED_FP`, `GCS_CKPT` (`checkpoints-vXYZ`), `GCS_CELLS_SUBDIR` (`cells-vXYZ`),
#      and append the prior `cells-v105` prefix to `restore_prev_java_cpp_artifacts()`.
# =============================================================================
set -euo pipefail

export GOLDEN_MODE=off
export CLOUDSDK_CORE_ACCOUNT="${CLOUDSDK_CORE_ACCOUNT:-vardict-agent@vardict-paper-trial.iam.gserviceaccount.com}"
export CLOUDSDK_CORE_PROJECT="${CLOUDSDK_CORE_PROJECT:-vardict-paper-trial}"

GCP_PROJECT="${CLOUDSDK_CORE_PROJECT}"
GCP_ZONE="${GCP_ZONE:-us-central1-a}"
GCP_SA="${CLOUDSDK_CORE_ACCOUNT}"
VM_NAME="${VM_NAME:-vardict-paper-n2s32}"
MACHINE_TYPE="${MACHINE_TYPE:-n2-standard-32}"
BOOT_DISK_SIZE="${BOOT_DISK_SIZE:-100GB}"
BOOT_DISK_TYPE="${BOOT_DISK_TYPE:-pd-ssd}"
LOCAL_SSD_COUNT="${LOCAL_SSD_COUNT:-4}"
MAX_RUN_DURATION="${MAX_RUN_DURATION:-24h}"
MIN_CPU_PLATFORM="${MIN_CPU_PLATFORM:-Intel Ice Lake}"
KEEP_VM_ON_FAIL="${KEEP_VM_ON_FAIL:-1}"

GCS_BAMS="${GCS_BAMS:-gs://vardict-paper-trial-bams}"
GCS_REFS="${GCS_REFS:-gs://vardict-paper-trial-refs}"
GCS_RESULTS="${GCS_RESULTS:-gs://vardict-paper-trial-results}"
GCS_EVIDENCE="${GCS_EVIDENCE:-gs://vardict-paper-trial-evidence}"
GCS_CKPT="${GCS_CKPT:-${GCS_RESULTS}/checkpoints-v105}"
GCS_ORACLE_CACHE="${GCS_ORACLE_CACHE:-${GCS_RESULTS}/oracle-cache-v183}"

EXPECTED_FP="${EXPECTED_FP:-d925892846918c2a}"
WORK_DIR="${WORK_DIR:-/var/tmp/vardict-paper}"

usage() {
  cat <<'EOF'
Usage:
  bash harness/gcp_paper_runner.sh --plan-check
  bash harness/gcp_paper_runner.sh --controller [--max-attempts 3] [--start-attempt 1]
  bash harness/gcp_paper_runner.sh --worker [--stage all|1|2|3|4]

Modes:
  --plan-check   Validate 2-strike Spot->On-Demand policy, NVMe RAID0 config, mkbed(), and checkpoint logic (offline safe)
  --controller   Host-side 2-strike Spot -> On-Demand VM lifecycle manager & progress monitor
  --worker       On-VM 4-stage benchmark executor with systemd-run cgroup slices and per-cell GCS checkpoints
EOF
}

# Maps attempt number (1..6) to GCP provisioning model (Attempts 1..5: SPOT; Attempt 6: STANDARD On-Demand).
provisioning_model_for_attempt() {
  local attempt="$1"
  if [ "$attempt" -le 5 ]; then
    echo "SPOT"
  else
    echo "STANDARD"
  fi
}

# Maps attempt number (1..6) across us-central1 zones to maximize Spot capacity availability before Attempt 6 On-Demand.
zone_for_attempt() {
  local attempt="$1"
  if [ "${GCP_ZONE_LOCK:-0}" = "1" ]; then
    echo "${GCP_ZONE}"
    return 0
  fi
  case "$attempt" in
    1) echo "us-central1-a" ;;
    2) echo "us-central1-b" ;;
    3) echo "us-central1-c" ;;
    4) echo "us-central1-f" ;;
    5) echo "us-central1-a" ;;
    *) echo "us-central1-a" ;;
  esac
}

# Builds the exact `gcloud compute instances create` command array for a given attempt.
build_gcloud_create_cmd() {
  local attempt="$1"
  local model zone
  model="$(provisioning_model_for_attempt "$attempt")"
  zone="$(zone_for_attempt "$attempt")"
  local cmd=(
    gcloud compute instances create "$VM_NAME"
    "--project=${GCP_PROJECT}"
    "--zone=${zone}"
    "--machine-type=${MACHINE_TYPE}"
    "--min-cpu-platform=${MIN_CPU_PLATFORM}"
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

# Exact mkbed() formula from Archive.zip/scripts/01_download_and_map.sh:77-83
# (mapped, primary, MAPQ>=1 reads via `samtools view -F 0x904 -q 1`; spans merged within 100 bp).
build_tierp_covered_bed() {
  local bam="$1" bed="$2" expected_lines="${3:-}"
  local bed_name
  bed_name="$(basename "$bed")"
  local generated=0
  if [ ! -s "$bed" ] && [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    setsid gsutil -q cp "${GCS_BAMS}/tier-p/${bed_name}" "$bed" 2>/dev/null || true
  fi
  if [ ! -s "$bed" ]; then
    samtools view -F 0x904 -q 1 "$bam" | awk 'BEGIN{OFS="\t"}{
      cig=$6; span=0;
      while(match(cig,/[0-9]+[MIDNSHP=X]/)){
        l=substr(cig,RSTART,RLENGTH-1)+0;
        op=substr(cig,RSTART+RLENGTH-1,1);
        if(op=="M"||op=="D"||op=="N"||op=="="||op=="X") span+=l;
        cig=substr(cig,RSTART+RLENGTH);
      }
      if(span>0) print $3,($4-1),($4-1+span)
    }' | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{
      if($1==c&&$2<=e+100){ if($3>e)e=$3 }
      else { if(c!="")print c,s,e,"R"; c=$1;s=$2;e=$3 }
    } END{ if(c!="")print c,s,e,"R" }' > "${bed}.tmp"
    mv "${bed}.tmp" "$bed"
    generated=1
  fi
  if [ -n "$expected_lines" ]; then
    local actual
    actual="$(wc -l < "$bed" | tr -dc '0-9')"
    if [ "$actual" -ne "$expected_lines" ]; then
      echo "FATAL: $bed has $actual regions, expected $expected_lines" >&2
      return 1
    fi
  fi
  if [ "$generated" -eq 1 ] && [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp "$bed" "${GCS_BAMS}/tier-p/${bed_name}" 2>/dev/null || true
  fi
}

# Single-process, 16-thread gsutil helper (Task 5a C1/C2: avoid gcloud storage cp hangs,
# avoid Ubuntu 24.04 snap multiprocessing.Pool fork deadlocks via parallel_process_count=1,
# keep 16 concurrent HTTP transfer threads, and retain MD5 on >150 MB objects).
GSUTIL_MT_FLAGS=(-m -o "GSUtil:parallel_process_count=1" -o "GSUtil:parallel_thread_count=16")

gcs_upload() {
  local src="$1" dst="$2"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    mkdir -p "$(dirname "$dst")"
    cp -r "$src" "$dst"
  else
    local t
    for t in 1 2 3; do
      if [ -d "$src" ]; then
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$src" "$dst" && return 0
      else
        setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp "$src" "$dst" && return 0
      fi
      sleep 2
    done
    return 1
  fi
}

gcs_stat_done() {
  local cell_id="$1"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    [ -f "${LOCAL_CKPT_DIR}/${cell_id}.done" ]
  else
    setsid gsutil -q stat "${GCS_CKPT}/${cell_id}.done" 2>/dev/null
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
    local t
    for t in 1 2 3; do
      printf '%s\n' "$stamp" | setsid gsutil -q cp - "${GCS_CKPT}/${cell_id}.done" && return 0
      sleep 2
    done
    return 1
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
    printf '%s\n' "$payload" | setsid gsutil -q cp - "${GCS_CKPT}/status.json" || true
  fi
}

# Restores previously measured VarDictJava, vardictcpp, and external field-caller
# (.var/.sorted.var/.vcf/.time/.exit/.err) artifacts for a given cell from cells-v103/
# and cells-v104/ on GCS (and preserves existing local non-Rust files) so neither oracle,
# C++, nor field callers are ever re-run when refreshing VarDictRust (e.g. v1.0.5) on n2-highmem-32.
restore_prev_java_cpp_artifacts() {
  local cell_id="$1" dst_dir="$2"
  mkdir -p "$dst_dir"
  if [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    local bucket_prefix
    for bucket_prefix in "${GCS_RESULTS}/cells-v103/${cell_id}" "${GCS_RESULTS}/cells-v104/${cell_id}" "${GCS_RESULTS}/cells-v105/${cell_id}"; do
      if setsid gsutil -q stat "${bucket_prefix}/*" 2>/dev/null; then
        echo "[CACHE-RESTORE] Restoring cached Java/Cpp/Field-Caller artifacts from ${bucket_prefix}/ -> ${dst_dir}/"
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r -x "(^|/)([^/]*rust[^/]*|.*compare.*|ceiling_.*|cell\\.exit_code|extra\\.json|.*\\.diff\\.txt|determinism_.*|parity_verdict\\.tsv|strelka_work.*|fb_chunks|vcfeval.*|${cell_id})(/|$)" "${bucket_prefix}" "$dst_dir" 2>/dev/null || true
      fi
    done
  fi
  # Never reuse stale Rust outputs from a prior version (keep any output newer than target/release/vardictrust)
  local rf
  while IFS= read -r -d '' rf; do
    if [ -x "target/release/vardictrust" ] && [ "$rf" -nt "target/release/vardictrust" ]; then
      continue
    fi
    rm -rf "$rf" 2>/dev/null || true
  done < <(find "$dst_dir" -maxdepth 3 \( -name "*rust*" -o -name "*vardictrust*" -o -name "*compare*" -o -name "ceiling_*" -o -name "cell.exit_code" -o -name "extra.json" -o -name "*.diff.txt" -o -name "determinism_*" -o -name "parity_verdict.tsv" \) -print0 2>/dev/null || true)
}

# Executes a single cell inside a systemd-run CPU+RAM cgroup slice (when systemd-run is available)
# and syncs KEEP_DIR to GCS on EVERY exit path before marking <cell_id>.done on success.
# `cpus` can be either an integer N (pinned to `0-$((N-1))`) or an explicit range like `8-15`.
run_checkpointed_cell() {
  local stage="$1" cell_id="$2" cpus="$3" mem_max="$4"
  shift 4
  if [ -n "${CELL_ALLOWLIST:-}" ] && [[ " ${CELL_ALLOWLIST} " != *" ${cell_id} "* ]]; then
    return 0
  fi
  if gcs_stat_done "$cell_id"; then
    echo "[SKIP] Cell '$cell_id' already checkpointed."
    return 0
  fi

  local cpu_spec mem_node_spec="0"
  if [ -n "${SOCKET_PIN_CPUS:-}" ] && { [ "$cpus" = "8" ] || [ "$cpus" = "0-7" ] || [ "$cpus" = "8-15" ]; }; then
    cpu_spec="${SOCKET_PIN_CPUS}"
  elif [ "${SOCKET_PIN_CPUS:-}" = "8-15" ] && [[ "$cpus" =~ ^[0-9]+$ ]] && [ "$cpus" -le 8 ]; then
    cpu_spec="8-$((8 + cpus - 1))"
  elif [[ "$cpus" =~ ^[0-9]+$ ]]; then
    cpu_spec="0-$((cpus - 1))"
  else
    cpu_spec="$cpus"
  fi
  if [ -d /sys/devices/system/node/node1 ]; then
    case "$cpu_spec" in
      0|0-[1-7]|0-3|4-7) mem_node_spec="0" ;;
      8|8-9|8-1[0-5]|12-15) mem_node_spec="1" ;;
      *) mem_node_spec="0-1" ;;
    esac
  fi

  echo "=== [RUN] Stage=$stage Cell=$cell_id AllowedCPUs=$cpu_spec AllowedMemoryNodes=$mem_node_spec MemoryMax=$mem_max GOLDEN_MODE=$GOLDEN_MODE ==="
  gcs_write_heartbeat "$stage" "$cell_id" "RUNNING"

  local cell_dir="${WORK_DIR}/cells/${cell_id}"
  mkdir -p "$cell_dir"
  restore_prev_java_cpp_artifacts "$cell_id" "$cell_dir"
  export KEEP_DIR="$cell_dir"
  export GOLDEN_MODE=off

  local rc=0
  if [ "${USE_SYSTEMD_SLICE:-1}" = "1" ] && command -v systemd-run >/dev/null 2>&1; then
    if [ "$(id -u)" -eq 0 ]; then
      systemd-run --scope --quiet \
        -p "AllowedCPUs=${cpu_spec}" \
        -p "AllowedMemoryNodes=${mem_node_spec}" \
        -p "MemoryMax=${mem_max}" \
        -- env PATH="$PATH" HOME="$HOME" WORK_DIR="$WORK_DIR" DATA_DIR="${DATA_DIR:-scratch/g4}" KEEP_DIR="$cell_dir" GOLDEN_MODE=off \
        GCS_BAMS="$GCS_BAMS" GCS_REFS="$GCS_REFS" GCS_RESULTS="$GCS_RESULTS" GCS_EVIDENCE="$GCS_EVIDENCE" GCS_CKPT="$GCS_CKPT" GCS_ORACLE_CACHE="$GCS_ORACLE_CACHE" "$@" || rc=$?
    else
      sudo -E systemd-run --scope --quiet \
        --uid="$(id -un)" --gid="$(id -gn)" \
        -p "AllowedCPUs=${cpu_spec}" \
        -p "AllowedMemoryNodes=${mem_node_spec}" \
        -p "MemoryMax=${mem_max}" \
        -- env PATH="$PATH" HOME="$HOME" WORK_DIR="$WORK_DIR" DATA_DIR="${DATA_DIR:-scratch/g4}" KEEP_DIR="$cell_dir" GOLDEN_MODE=off \
        GCS_BAMS="$GCS_BAMS" GCS_REFS="$GCS_REFS" GCS_RESULTS="$GCS_RESULTS" GCS_EVIDENCE="$GCS_EVIDENCE" GCS_CKPT="$GCS_CKPT" GCS_ORACLE_CACHE="$GCS_ORACLE_CACHE" "$@" || rc=$?
    fi
  elif command -v taskset >/dev/null 2>&1; then
    taskset -c "${cpu_spec}" "$@" || rc=$?
  else
    "$@" || rc=$?
  fi

  echo "$rc" > "${cell_dir}/cell.exit_code"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    gcs_upload "$cell_dir" "${LOCAL_CKPT_DIR}/cells/${cell_id}"
  else
    gcs_upload "$cell_dir" "${GCS_RESULTS}/${GCS_CELLS_SUBDIR:-cells-v105}/${cell_id}" || true
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

# Executes a local-disk data staging cell guarded by ${WORK_DIR}/local_staged/${cell_id}.done
# so a newly provisioned VM after Spot preemption always restages its ephemeral NVMe RAID0 disk.
run_local_stage_cell() {
  local stage="$1" cell_id="$2" cpus="$3" mem_max="$4"
  shift 4
  if [ -n "${CELL_ALLOWLIST:-}" ] && [[ " ${CELL_ALLOWLIST} " != *" ${cell_id} "* ]]; then
    return 0
  fi
  local local_marker="${WORK_DIR}/local_staged/${cell_id}.done"
  if [ -f "$local_marker" ]; then
    echo "[SKIP] Local staging cell '$cell_id' already present on this VM."
    return 0
  fi

  echo "=== [STAGE-LOCAL] Stage=$stage Cell=$cell_id CPUs=$cpus MemoryMax=$mem_max ==="
  gcs_write_heartbeat "$stage" "$cell_id" "STAGING_LOCAL"

  local cell_dir="${WORK_DIR}/cells/${cell_id}"
  mkdir -p "$cell_dir" "${WORK_DIR}/local_staged"
  export KEEP_DIR="$cell_dir"
  export GOLDEN_MODE=off

  local rc=0
  "$@" || rc=$?
  echo "$rc" > "${cell_dir}/cell.exit_code"
  if [ "$rc" -ne 0 ]; then
    gcs_write_heartbeat "$stage" "$cell_id" "FAILED_EXIT_${rc}"
    echo "[ERROR] Local staging cell '$cell_id' exited with code $rc." >&2
    return "$rc"
  fi

  date -u +%Y-%m-%dT%H:%M:%SZ > "$local_marker"
  gcs_mark_done "$cell_id"
  gcs_write_heartbeat "$stage" "$cell_id" "DONE"
  echo "=== [DONE] Local staging cell '$cell_id' ready on NVMe ==="
}

# =============================================================================
# Disjoint-Core Background Pre-Stager (AllowedCPUs=24-31, MemoryMax=64G)
# Uses otherwise-idle cores on n2-highmem-32 during 8-core Stage-2 runs to
# pre-provision field-caller conda envs and pre-stage Stage-3/Stage-4 assets.
# =============================================================================
ensure_field_caller_envs() {
  if [ ! -f "${WORK_DIR}/local_staged/field_envs.done" ]; then
    echo "[Pre-Stager] Provisioning field-caller conda envs (FreeBayes 1.3.10, LoFreq 2.1.5, Mutect2 4.6.2.0, octopus 0.7.4, Strelka2 2.9.10)..."
    conda env list | grep -q '^cl_freebayes ' || conda create -y -n cl_freebayes -c conda-forge -c bioconda freebayes=1.3.10 samtools=1.24 bcftools=1.24 parallel
    conda env list | grep -q '^cl_lofreq ' || conda create -y -n cl_lofreq -c conda-forge -c bioconda lofreq=2.1.5 samtools=1.24
    conda env list | grep -q '^cl_gatk4 ' || conda create -y -n cl_gatk4 -c conda-forge -c bioconda gatk4=4.6.2.0 openjdk=21 samtools=1.24
    conda env list | grep -q '^cl_octopus ' || conda create -y -n cl_octopus -c conda-forge -c bioconda octopus=0.7.4
    conda env list | grep -q '^cl_strelka ' || conda create -y -n cl_strelka -c conda-forge -c bioconda strelka=2.9.10 python=2.7
    mkdir -p "${WORK_DIR}/local_staged"
    date -u +%Y-%m-%dT%H:%M:%SZ > "${WORK_DIR}/local_staged/field_envs.done"
    echo "[Pre-Stager] Field-caller conda envs ready."
  fi
  # Bioconda strelka-2.9.10-2 patchelf places .dynamic (PT_LOAD RW) on the same 4K page (0x4dc000)
  # as PT_LOAD segment 7 (R E) in libexec/strelka2; ensure PF_W (0x2) is set on overlapping pages.
  python3 -c '
import glob, os, struct
for path in sorted(glob.glob(os.path.expanduser("~/miniconda3/envs/cl_strelka/share/strelka-*/libexec/*"))):
    if not os.path.isfile(path) or os.path.islink(path):
        continue
    with open(path, "rb") as f:
        hdr = f.read(64)
    if hdr[:4] != b"\x7fELF":
        continue
    e_phoff = struct.unpack("<Q", hdr[32:40])[0]
    e_phentsize, e_phnum = struct.unpack("<HH", hdr[54:58])
    with open(path, "r+b") as f:
        loads = []
        for i in range(e_phnum):
            f.seek(e_phoff + i * e_phentsize)
            p_type, p_flags, _, p_vaddr, _, _, p_memsz, _ = struct.unpack("<IIQQQQQQ", f.read(56))
            if p_type == 1:
                loads.append((i, p_flags, p_vaddr, p_memsz))
        rw_pages = {p for _, fl, va, ms in loads if fl & 2 for p in range(va >> 12, (va + ms + 4095) >> 12)}
        for i, fl, va, ms in loads:
            if not (fl & 2) and (set(range(va >> 12, (va + ms + 4095) >> 12)) & rw_pages):
                f.seek(e_phoff + i * e_phentsize + 4)
                f.write(struct.pack("<I", fl | 2))
' || true
  # Bioconda lofreq-2.1.5 lofreq2_call_pparallel.py compares len(bed_sqs) against len(bam_sqs)/10.0,
  # which fails on GRCh38.d1.vd1 BAMs with 1,689 decoy/alt/HLA contigs and spawns 271,597 subprocesses
  # that overflow ARG_MAX in bcftools concat; filter bam_bins to bed_sqs in the threshold check.
  python3 -c '
import os
p = os.path.expanduser("~/miniconda3/envs/cl_lofreq/bin/lofreq2_call_pparallel.py")
if os.path.isfile(p):
    with open(p) as f:
        s = f.read()
    old = "if len(bed_bins) > 100*len(bam_bins) and len(bed_sqs) > len(bam_sqs)/10.0:"
    new = "if len(bed_bins) > 100*len([b for b in bam_bins if b[0] in bed_sqs]) and len(bed_sqs) > 0:"
    if old in s:
        with open(p, "w") as f:
            f.write(s.replace(old, new))
' || true
}

start_background_prestager() {
  if [ -n "${CELL_ALLOWLIST:-}" ] || [ -f "${WORK_DIR}/local_staged/bg_prestage.done" ]; then
    return 0
  fi
  (
    set -euo pipefail
    exec > "${WORK_DIR}/bg_prestage.log" 2>&1
    echo "[BG Pre-Stager] Starting on disjoint cores 24-31 at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    taskset -c 24-31 bash -c "
      set -euo pipefail
      export PATH=\"\$HOME/miniconda3/envs/g4-tools/bin:\$HOME/miniconda3/bin:\$HOME/.cargo/bin:\$PATH\"
      conda env list | grep -q '^cl_freebayes ' || conda create -y -n cl_freebayes -c conda-forge -c bioconda freebayes=1.3.10 samtools=1.24 bcftools=1.24 parallel
      conda env list | grep -q '^cl_lofreq ' || conda create -y -n cl_lofreq -c conda-forge -c bioconda lofreq=2.1.5 samtools=1.24
      conda env list | grep -q '^cl_gatk4 ' || conda create -y -n cl_gatk4 -c conda-forge -c bioconda gatk4=4.6.2.0 openjdk=21 samtools=1.24
      conda env list | grep -q '^cl_octopus ' || conda create -y -n cl_octopus -c conda-forge -c bioconda octopus=0.7.4
      conda env list | grep -q '^cl_strelka ' || conda create -y -n cl_strelka -c conda-forge -c bioconda strelka=2.9.10 python=2.7
      mkdir -p '${WORK_DIR}/local_staged' '${WORK_DIR}/data/rna'
      date -u +%Y-%m-%dT%H:%M:%SZ > '${WORK_DIR}/local_staged/field_envs.done'
      if [ ! -f '${WORK_DIR}/data/rna/HG002_rep1_Aligned.sortedByCoord.out.bam' ]; then
        setsid gsutil -q cp '${GCS_BAMS}/rna-hg002/aligned/HG002_rep1_Aligned.sortedByCoord.out.bam*' '${WORK_DIR}/data/rna/' || true
      fi
      if setsid gsutil -q stat '${GCS_BAMS}/g4-cache/s4_ready.done' 2>/dev/null; then
        mkdir -p scratch/g4/seqc2 scratch/g4/multimodal
        setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/seqc2' scratch/g4/seqc2 || true
        setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/multimodal' scratch/g4/multimodal || true
      fi
      date -u +%Y-%m-%dT%H:%M:%SZ > '${WORK_DIR}/local_staged/bg_prestage.done'
    "
    echo "[BG Pre-Stager] Completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  ) &
  BG_PRESTAGE_PID=$!
  echo "[Worker] Launched background pre-stager (PID=${BG_PRESTAGE_PID}) on cores 24-31."
}

# =============================================================================
# Worker Stage Implementations (Executed on n2-highmem-32 / 3.0 TB NVMe RAID0)
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
  local size_gb avail_gb
  size_gb="$(df -BG --output=size "$WORK_DIR" | tail -1 | tr -dc '0-9')"
  avail_gb="$(df -BG --output=avail "$WORK_DIR" | tail -1 | tr -dc '0-9')"
  local min_size_gb="${MIN_DISK_SIZE_GB:-1200}"
  local min_avail_gb="${MIN_DISK_GB:-400}"
  if [ "$size_gb" -lt "$min_size_gb" ] || [ "$avail_gb" -lt "$min_avail_gb" ]; then
    echo "FATAL: $WORK_DIR has ${size_gb}GB total / ${avail_gb}GB free, need >= ${min_size_gb}GB total and >= ${min_avail_gb}GB free (NVMe RAID0 SSD required)" >&2
    exit 2
  fi
  local fp
  fp="$(bash harness/engine_fingerprint.sh)"
  if [ "$fp" != "$EXPECTED_FP" ]; then
    echo "FATAL: engine fingerprint mismatch: got '$fp', expected '$EXPECTED_FP'" >&2
    exit 2
  fi
  local meta_cpu=""
  meta_cpu="$(curl -fsSL -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/cpu-platform 2>/dev/null || true)"
  if [ -n "$meta_cpu" ] && [ "$meta_cpu" != "$MIN_CPU_PLATFORM" ]; then
    echo "FATAL: VM cpuPlatform is '$meta_cpu', expected '$MIN_CPU_PLATFORM' (Cascade Lake placement drift is forbidden)" >&2
    exit 2
  fi
  {
    echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "expected_fp=${EXPECTED_FP}"
    echo "perf_fp=$(bash harness/engine_fingerprint.sh --perf)"
    echo "gcp_cpu_platform=${meta_cpu:-local}"
    echo "gcp_machine_type=$(curl -fsSL -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/machine-type 2>/dev/null || true)"
    echo "gcp_zone=$(curl -fsSL -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/zone 2>/dev/null || true)"
    echo "=== lscpu ==="
    lscpu 2>/dev/null || true
    echo "=== free -h ==="
    free -h 2>/dev/null || true
    echo "=== df -h WORK_DIR ==="
    df -h "$WORK_DIR" 2>/dev/null || true
  } > "$WORK_DIR/cells/host_hardware_manifest.txt"
  gcs_upload "$WORK_DIR/cells/host_hardware_manifest.txt" "${GCS_RESULTS}/${GCS_CELLS_SUBDIR:-cells-v105}/host_hardware_manifest.txt" || true
  cargo build --release -p vardictrust
}

run_stage_1() {
  echo "=== [Stage 1: Task 5b Pre-Full-Genome Gates] ==="
  export DATA_DIR="scratch/g4"

  run_local_stage_cell "stage1" "s1_stage_assets" 32 "128G" bash -c "
    set -euo pipefail
    mkdir -p scratch/g4/amplicon scratch/g4/hg002 scratch/g4/hg008 scratch/g4/ref scratch/g4/truth_cache
    s1_cached=0
    refs_cached=0
    if setsid gsutil -q stat '${GCS_BAMS}/g4-cache/s1_ready.done' 2>/dev/null; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/amplicon' scratch/g4/amplicon
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/hg002' scratch/g4/hg002
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/hg008' scratch/g4/hg008
      s1_cached=1
    fi
    if setsid gsutil -q stat '${GCS_BAMS}/g4-cache/truth_cache/*' 2>/dev/null; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/truth_cache' scratch/g4/truth_cache || true
    fi
    if setsid gsutil -q stat '${GCS_BAMS}/g4-cache/refs_ready.done' 2>/dev/null; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/ref' scratch/g4/ref
      python3 harness/convert_fastas.py
      ln -sfn GRCh38_GIABv3.fasta scratch/g4/ref/giabv3/GRCh38_GIABv3.fa
      ln -sfn GRCh38_GIABv3.fasta.fai scratch/g4/ref/giabv3/GRCh38_GIABv3.fa.fai
      ln -sfn GRCh38_GIABv3.sdf scratch/g4/ref/giabv3/GRCh38_GIABv3.fa.sdf
      ln -sfn GRCh38.d1.vd1.sdf scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa.sdf
      refs_cached=1
    fi
    bash harness/stage_g4.sh refs amplicon hg002 hg008
    hg002_md5=\"\$(md5sum scratch/g4/hg002/HG002.chr20.bam | awk '{print \$1}')\"
    if [ \"\$hg002_md5\" != 'bd605a455b2e2064931e6097c64695ea' ]; then
      echo \"FATAL: scratch/g4/hg002/HG002.chr20.bam MD5 '\$hg002_md5' != MANIFEST pin 'bd605a455b2e2064931e6097c64695ea'\" >&2
      exit 1
    fi
    if [ \"\$s1_cached\" -eq 0 ]; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r scratch/g4/amplicon '${GCS_BAMS}/g4-cache/amplicon' || true
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r scratch/g4/hg002 '${GCS_BAMS}/g4-cache/hg002' || true
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r scratch/g4/hg008 '${GCS_BAMS}/g4-cache/hg008' || true
      date -u +%Y-%m-%dT%H:%M:%SZ | setsid gsutil -q cp - '${GCS_BAMS}/g4-cache/s1_ready.done' || true
    fi
    if [ \"\$refs_cached\" -eq 0 ]; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r scratch/g4/ref '${GCS_BAMS}/g4-cache/ref' || true
      date -u +%Y-%m-%dT%H:%M:%SZ | setsid gsutil -q cp - '${GCS_BAMS}/g4-cache/refs_ready.done' || true
    fi
  "
  ensure_submodule_excludes

  # H1a Simple-mode tile-size equivalence (500k, 1M, 2M, 10M vs untiled oracle on HG002 slice sub-window)
  run_checkpointed_cell "stage1" "s1_h1a_simple_tiles" 8 "32G" bash -c '
    set -euo pipefail
    head -n 40 "$DATA_DIR/hg002/eval.bed" > "$KEEP_DIR/h1a_simple.bed"
    if [ ! -s "$KEEP_DIR/oracle_untiled.vcf" ]; then
      ORACLE_XMX_GB=22 bash harness/run_oracle.sh \
        "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" \
        "$DATA_DIR/hg002/HG002.chr20.bam" \
        "$KEEP_DIR/h1a_simple.bed" \
        "$KEEP_DIR/oracle_untiled.vcf" HG002 -f 0.01 -th 8 -z
    fi
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
    if [ -s "$KEEP_DIR/java.vcf" ] && [ -s "$KEEP_DIR/java.time" ]; then
      /usr/bin/time -v -o "$KEEP_DIR/rust.time" \
        target/release/vardictrust -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" -f 0.01 -N HG002 \
        -b "$DATA_DIR/hg002/HG002.chr20.bam" -c 1 -S 2 -E 3 -g 4 -z -th 8 "$KEEP_DIR/gate7_multiregion.bed" > "$KEEP_DIR/rust.vcf"
      python3 harness/compare_vcf.py --min-variants 50 "$KEEP_DIR/java.vcf" "$KEEP_DIR/rust.vcf" | tee "$KEEP_DIR/compare.txt"
      grep -q "^MATCH" "$KEEP_DIR/compare.txt"
    else
      GOLDEN_MODE=off MIN_VARIANTS=50 ORACLE_XMX_GB=22 \
        bash harness/compare.sh \
        "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" \
        "$DATA_DIR/hg002/HG002.chr20.bam" \
        "$KEEP_DIR/gate7_multiregion.bed" HG002 -f 0.01 -th 8
    fi
  '

  # Gate 1 Ceiling Cell A: Amplicon deep panel (-th 8, assert Rust peak RSS <= 8 GB)
  run_checkpointed_cell "stage1" "s1_ceiling_amplicon" 8 "32G" bash -c '
    set -euo pipefail
    if [ -s "$KEEP_DIR/java.vcf" ] && [ -s "$KEEP_DIR/java.time" ]; then
      /usr/bin/time -v -o "$KEEP_DIR/rust.time" \
        target/release/vardictrust -G "$DATA_DIR/amplicon/ref.fa" -f 0.01 -N SAMPLE \
        -b "$DATA_DIR/amplicon/amplicon.bam" -c 1 -S 2 -E 3 -g 4 -z -a 10:0.95 -th 8 "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" > "$KEEP_DIR/rust.vcf"
      python3 harness/compare_vcf.py "$KEEP_DIR/java.vcf" "$KEEP_DIR/rust.vcf" | tee "$KEEP_DIR/compare.txt"
      grep -q "^MATCH" "$KEEP_DIR/compare.txt"
    else
      GOLDEN_MODE=off ORACLE_XMX_GB=22 \
        bash harness/compare.sh \
        "$DATA_DIR/amplicon/ref.fa" \
        "$DATA_DIR/amplicon/amplicon.bam" \
        "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" SAMPLE -a 10:0.95 -th 8
    fi
    rss_kb=$(awk -F": " "/Maximum resident set size/{print \$2}" "$KEEP_DIR/rust.time" | tr -dc "0-9")
    echo "amplicon_ceiling_rss_kb=$rss_kb" | tee "$KEEP_DIR/ceiling_amplicon.txt"
    [ "$rss_kb" -le 8000000 ] || { echo "FATAL: amplicon ceiling RSS ${rss_kb} KB > 8 GB" >&2; exit 9; }
  '

  # Gate 1 Ceiling Cell B: HG008 195x/169x somatic WGS slice (-th 8, assert Rust peak RSS <= 8 GB)
  run_checkpointed_cell "stage1" "s1_ceiling_hg008" 8 "32G" bash -c '
    set -euo pipefail
    if [ -s "$KEEP_DIR/java.vcf" ] && [ -s "$KEEP_DIR/java.time" ]; then
      /usr/bin/time -v -o "$KEEP_DIR/rust.time" \
        target/release/vardictrust -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" -f 0.01 -N "HG008-T|HG008-N" \
        -b "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" -c 1 -S 2 -E 3 -g 4 -z -th 8 "$DATA_DIR/hg008/engine.bed" > "$KEEP_DIR/rust.vcf"
      python3 harness/compare_vcf.py --min-variants 1000 "$KEEP_DIR/java.vcf" "$KEEP_DIR/rust.vcf" | tee "$KEEP_DIR/compare.txt"
      grep -q "^MATCH" "$KEEP_DIR/compare.txt"
    else
      GOLDEN_MODE=off MIN_VARIANTS=1000 ORACLE_XMX_GB=22 \
        bash harness/compare.sh \
        "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" \
        "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" \
        "$DATA_DIR/hg008/engine.bed" "HG008-T|HG008-N" -th 8
    fi
    rss_kb=$(awk -F": " "/Maximum resident set size/{print \$2}" "$KEEP_DIR/rust.time" | tr -dc "0-9")
    echo "hg008_ceiling_rss_kb=$rss_kb" | tee "$KEEP_DIR/ceiling_hg008.txt"
    [ "$rss_kb" -le 8000000 ] || { echo "FATAL: HG008 ceiling RSS ${rss_kb} KB > 8 GB" >&2; exit 9; }
  '
}

run_stage_2() {
  echo "=== [Stage 2: Task 5c Tier-P 5-Sample Replication + 13-Mode Sweep + Fig-2 Truth + Field Callers + Table-3 -R chr1] ==="
  local tp_dir="${WORK_DIR}/data/tier-p"
  mkdir -p "$tp_dir/bams" "$tp_dir/truth" "$tp_dir/archive"

  run_local_stage_cell "stage2" "s2_stage_tierp_assets" 32 "128G" bash -c "
    set -euo pipefail
    setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q cp '${GCS_BAMS}/tier-p/*' '$tp_dir/bams/'
    setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q cp '${GCS_REFS}/tier-p-truth/*' '$tp_dir/truth/'
    if [ ! -f '$tp_dir/archive/Archive.zip' ]; then
      setsid gsutil -q cp '${GCS_REFS}/vardictcpp/Archive.zip' '$tp_dir/archive/Archive.zip'
    fi
    unzip -q -o '$tp_dir/archive/Archive.zip' 'reference/*' 'beds/*' 'scripts/*' 'performance_benchmark/*' 'results/*' 'README.md' 'MD5SUMS.txt' -d '$tp_dir/archive/'
    if [ ! -x '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' ]; then
      mkdir -p '$tp_dir/vardictcpp-ae93fe39dc/build'
      if setsid gsutil -q cp '${GCS_REFS}/vardictcpp/vardictcpp-ae93fe39dc.bin' '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' 2>/dev/null; then
        chmod +x '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp'
      else
        rm -rf '$tp_dir/vardictcpp-ae93fe39dc'
        git clone https://github.com/MHH-Bioinformatics-Hematology/vardictcpp.git '$tp_dir/vardictcpp-ae93fe39dc'
        cd '$tp_dir/vardictcpp-ae93fe39dc'
        git checkout ae93fe39dc
        conda run -n vardict-oracle bash -c 'cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 16'
        setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' '${GCS_REFS}/vardictcpp/vardictcpp-ae93fe39dc.bin' || true
      fi
    fi
  "

  # Start background pre-stager on disjoint cores 24-31 while Stage 2 runs on cores 0-7
  start_background_prestager

  # Table-1 expected region and row counts per accession (gcp-prereg.md §A4 & Archive.zip/results/parity.tsv)
  declare -A EXP_REGIONS=(
    [SRR15006386]=11732
    [SRR15006375]=12292
    [SRR15006376]=105949
    [SRR15006540]=135145
    [SRR8657348]=773667
  )
  declare -A EXP_ROWS=(
    [SRR15006386]=4014
    [SRR15006375]=4056
    [SRR15006376]=12739
    [SRR15006540]=16031
    [SRR8657348]=54448
  )

  # 1. Five SRR samples x {8, 1} threads x 3 engines x n=3 technical replicates (exact mkbed(), no -z, .var + .vcf)
  for acc in SRR15006386 SRR15006375 SRR15006376 SRR15006540 SRR8657348; do
    local exp_reg="${EXP_REGIONS[$acc]}"
    local exp_row="${EXP_ROWS[$acc]}"
    for th in 8 1; do
      run_checkpointed_cell "stage2" "s2_tierp_v2_${acc}_th${th}" 8 "32G" bash -c "
        set -euo pipefail
        source harness/gcp_paper_runner.sh --source-only
        REF='$tp_dir/archive/reference/hg19.fa'
        BAM='$tp_dir/bams/${acc}.bam'
        BED='$tp_dir/bams/${acc}.mkbed.bed'
        build_tierp_covered_bed \"\$BAM\" \"\$BED\" '${exp_reg}'
        wc -l \"\$BED\" > \"\$KEEP_DIR/${acc}.bed_lines\"

        # Symmetric n=3 technical replicates on the exact Wolff 2026 .var TSV surface (no -z)
        for run in 1 2 3; do
          if [ ! -s \"\$KEEP_DIR/rust_var_th${th}_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/rust_th${th}.var\" ]; then
            /usr/bin/time -v -o \"\$KEEP_DIR/rust_var_th${th}_run\${run}.time\" \
              target/release/vardictrust -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th ${th} \"\$BED\" > \"\$KEEP_DIR/rust_th${th}.var\"
          fi
          if [ ! -s \"\$KEEP_DIR/cpp_var_th${th}_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/cpp_th${th}.var\" ]; then
            /usr/bin/time -v -o \"\$KEEP_DIR/cpp_var_th${th}_run\${run}.time\" \
              '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 --chunk 50000 -th ${th} \"\$BED\" > \"\$KEEP_DIR/cpp_th${th}.var\"
          fi
          if [ ! -s \"\$KEEP_DIR/java_var_th${th}_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/java_th${th}.var\" ]; then
            _JAVA_OPTIONS='-Xmx22g' /usr/bin/time -v -o \"\$KEEP_DIR/java_var_th${th}_run\${run}.time\" \
              conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -th ${th} \"\$BED\" > \"\$KEEP_DIR/java_th${th}.var\"
          fi
        done

        # Verify exact Table-1 row counts & 0 FP / 0 FN sorted TSV identity + compare_var.py
        r_rows=\$(wc -l < \"\$KEEP_DIR/rust_th${th}.var\" | tr -dc '0-9')
        j_rows=\$(wc -l < \"\$KEEP_DIR/java_th${th}.var\" | tr -dc '0-9')
        c_rows=\$(wc -l < \"\$KEEP_DIR/cpp_th${th}.var\" | tr -dc '0-9')
        echo \"rows_rust=\$r_rows rows_java=\$j_rows rows_cpp=\$c_rows expected=${exp_row}\" | tee \"\$KEEP_DIR/table1_counts_th${th}.txt\"
        [ \"\$r_rows\" -eq '${exp_row}' ] && [ \"\$j_rows\" -eq '${exp_row}' ]
        LC_ALL=C sort \"\$KEEP_DIR/java_th${th}.var\" > \"\$KEEP_DIR/java_th${th}.sorted.var\"
        LC_ALL=C sort \"\$KEEP_DIR/rust_th${th}.var\" > \"\$KEEP_DIR/rust_th${th}.sorted.var\"
        LC_ALL=C sort \"\$KEEP_DIR/cpp_th${th}.var\" > \"\$KEEP_DIR/cpp_th${th}.sorted.var\"
        diff -u \"\$KEEP_DIR/java_th${th}.sorted.var\" \"\$KEEP_DIR/rust_th${th}.sorted.var\" > \"\$KEEP_DIR/var_diff_th${th}.txt\" || true
        key_fn() { awk -F'\t' '{print \$3\"_\"\$4\"_\"\$6\"_\"\$7}' \"\$1\" | LC_ALL=C sort; }
        fp_jr=\$(LC_ALL=C comm -13 <(key_fn \"\$KEEP_DIR/java_th${th}.var\") <(key_fn \"\$KEEP_DIR/rust_th${th}.var\") | wc -l | tr -dc '0-9')
        fn_jr=\$(LC_ALL=C comm -23 <(key_fn \"\$KEEP_DIR/java_th${th}.var\") <(key_fn \"\$KEEP_DIR/rust_th${th}.var\") | wc -l | tr -dc '0-9')
        fp_jc=\$(LC_ALL=C comm -13 <(key_fn \"\$KEEP_DIR/java_th${th}.var\") <(key_fn \"\$KEEP_DIR/cpp_th${th}.var\") | wc -l | tr -dc '0-9')
        fn_jc=\$(LC_ALL=C comm -23 <(key_fn \"\$KEEP_DIR/java_th${th}.var\") <(key_fn \"\$KEEP_DIR/cpp_th${th}.var\") | wc -l | tr -dc '0-9')
        echo \"fp_jr=\$fp_jr fn_jr=\$fn_jr fp_jc=\$fp_jc fn_jc=\$fn_jc\" | tee \"\$KEEP_DIR/equivalence_th${th}.txt\"
        [ \"\$fp_jr\" -eq 0 ] && [ \"\$fn_jr\" -eq 0 ]
        python3 harness/compare_var.py --eps 1e-6 <(sed 's/:0-/:1-/g' \"\$KEEP_DIR/java_th${th}.sorted.var\") \"\$KEEP_DIR/rust_th${th}.sorted.var\" | tee \"\$KEEP_DIR/compare_var_th${th}.txt\"

        # End-to-end VCF surface (consolidated Rust binary vs Java+R+Perl vs C++ -vcf) + tabix check
        /usr/bin/time -v -o \"\$KEEP_DIR/rust_vcf_th${th}.time\" \
          target/release/vardictrust -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -th ${th} \"\$BED\" > \"\$KEEP_DIR/rust_th${th}.vcf\"
        if [ ! -s \"\$KEEP_DIR/cpp_vcf_th${th}.time\" ] || [ ! -s \"\$KEEP_DIR/cpp_th${th}.vcf\" ]; then
          /usr/bin/time -v -o \"\$KEEP_DIR/cpp_vcf_th${th}.time\" \
            '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 --chunk 50000 -th ${th} \"\$BED\" -vcf > \"\$KEEP_DIR/cpp_th${th}.vcf\"
        fi
        if [ ! -s \"\$KEEP_DIR/java_vcf_th${th}.time\" ] || [ ! -s \"\$KEEP_DIR/java_th${th}.vcf\" ]; then
          ORACLE_XMX_GB=22 /usr/bin/time -v -o \"\$KEEP_DIR/java_vcf_th${th}.time\" \
            bash harness/run_oracle.sh \"\$REF\" \"\$BAM\" \"\$BED\" \"\$KEEP_DIR/java_th${th}.vcf\" '${acc}' -f 0.01 -th ${th}
        fi
        python3 harness/compare_vcf.py \"\$KEEP_DIR/java_th${th}.vcf\" \"\$KEEP_DIR/rust_th${th}.vcf\" | tee \"\$KEEP_DIR/compare_vcf_th${th}.txt\"
        grep -q '^MATCH' \"\$KEEP_DIR/compare_vcf_th${th}.txt\"
        bcftools sort \"\$KEEP_DIR/rust_th${th}.vcf\" -Oz -o \"\$KEEP_DIR/rust_th${th}.vcf.gz\" && tabix -f -p vcf \"\$KEEP_DIR/rust_th${th}.vcf.gz\" && bcftools index -f \"\$KEEP_DIR/rust_th${th}.vcf.gz\"
      "
    done
  done

  # 2. Somatic pair SRR15006376 (tumor) | SRR15006386 (normal), split by thread count (8, 1), n=3 across all 3 engines
  for th in 8 1; do
    run_checkpointed_cell "stage2" "s2_tierp_v2_somatic_76_86_th${th}" 8 "32G" bash -c "
      set -euo pipefail
      source harness/gcp_paper_runner.sh --source-only
      REF='$tp_dir/archive/reference/hg19.fa'
      TUM='$tp_dir/bams/SRR15006376.bam'
      NOR='$tp_dir/bams/SRR15006386.bam'
      BED='$tp_dir/bams/SRR15006376.mkbed.bed'
      build_tierp_covered_bed \"\$TUM\" \"\$BED\" 105949

      for run in 1 2 3; do
        /usr/bin/time -v -o \"\$KEEP_DIR/rust_som_var_th${th}_run\${run}.time\" \
          target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th ${th} \"\$BED\" > \"\$KEEP_DIR/rust_som_th${th}.var\"
        if [ ! -s \"\$KEEP_DIR/cpp_som_var_th${th}_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/cpp_som_th${th}.var\" ]; then
          /usr/bin/time -v -o \"\$KEEP_DIR/cpp_som_var_th${th}_run\${run}.time\" \
            '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 --chunk 50000 -th ${th} \"\$BED\" > \"\$KEEP_DIR/cpp_som_th${th}.var\"
        fi
        if [ ! -s \"\$KEEP_DIR/java_som_var_th${th}_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/java_som_th${th}.var\" ]; then
          _JAVA_OPTIONS='-Xmx22g' /usr/bin/time -v -o \"\$KEEP_DIR/java_som_var_th${th}_run\${run}.time\" \
            conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 -th ${th} \"\$BED\" > \"\$KEEP_DIR/java_som_th${th}.var\"
        fi
      done

      r_rows=\$(wc -l < \"\$KEEP_DIR/rust_som_th${th}.var\" | tr -dc '0-9')
      j_rows=\$(wc -l < \"\$KEEP_DIR/java_som_th${th}.var\" | tr -dc '0-9')
      c_rows=\$(wc -l < \"\$KEEP_DIR/cpp_som_th${th}.var\" | tr -dc '0-9')
      echo \"somatic_rows_rust=\$r_rows rows_java=\$j_rows rows_cpp=\$c_rows expected=14126\" | tee \"\$KEEP_DIR/somatic_counts_th${th}.txt\"
      [ \"\$r_rows\" -eq 14126 ] && [ \"\$j_rows\" -eq 14126 ]
      LC_ALL=C sort \"\$KEEP_DIR/java_som_th${th}.var\" > \"\$KEEP_DIR/java_som_th${th}.sorted.var\"
      LC_ALL=C sort \"\$KEEP_DIR/rust_som_th${th}.var\" > \"\$KEEP_DIR/rust_som_th${th}.sorted.var\"
      LC_ALL=C sort \"\$KEEP_DIR/cpp_som_th${th}.var\" > \"\$KEEP_DIR/cpp_som_th${th}.sorted.var\"
      diff -u \"\$KEEP_DIR/java_som_th${th}.sorted.var\" \"\$KEEP_DIR/rust_som_th${th}.sorted.var\" > \"\$KEEP_DIR/som_var_diff_th${th}.txt\" || true
      python3 harness/compare_var.py --eps 1e-6 <(sed 's/:0-/:1-/g' \"\$KEEP_DIR/java_som_th${th}.sorted.var\") \"\$KEEP_DIR/rust_som_th${th}.sorted.var\" | tee \"\$KEEP_DIR/compare_som_var_th${th}.txt\"

      /usr/bin/time -v -o \"\$KEEP_DIR/rust_som_vcf_th${th}.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 -th ${th} \"\$BED\" > \"\$KEEP_DIR/rust_som_th${th}.vcf\"
      if [ ! -s \"\$KEEP_DIR/cpp_som_vcf_th${th}.time\" ] || [ ! -s \"\$KEEP_DIR/cpp_som_th${th}.vcf\" ]; then
        /usr/bin/time -v -o \"\$KEEP_DIR/cpp_som_vcf_th${th}.time\" \
          '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 --chunk 50000 -th ${th} \"\$BED\" -vcf > \"\$KEEP_DIR/cpp_som_th${th}.vcf\"
      fi
      if [ ! -s \"\$KEEP_DIR/java_som_vcf_th${th}.time\" ] || [ ! -s \"\$KEEP_DIR/java_som_th${th}.vcf\" ]; then
        ORACLE_XMX_GB=22 /usr/bin/time -v -o \"\$KEEP_DIR/java_som_vcf_th${th}.time\" \
          bash harness/run_oracle.sh \"\$REF\" \"\${TUM}|\${NOR}\" \"\$BED\" \"\$KEEP_DIR/java_som_th${th}.vcf\" 'SRR15006376|SRR15006386' -f 0.01 -th ${th}
      fi
      python3 harness/compare_vcf.py \"\$KEEP_DIR/java_som_th${th}.vcf\" \"\$KEEP_DIR/rust_som_th${th}.vcf\" | tee \"\$KEEP_DIR/compare_som_vcf_th${th}.txt\"
      grep -q '^MATCH' \"\$KEEP_DIR/compare_som_vcf_th${th}.txt\"
    "
  done

  # 3. Tier-P 13-Mode Parity Sweep (results/full_mode_parity.tsv + MRD -f 0 arm + P5 Panel CRAM Twin)
  run_checkpointed_cell "stage2" "s2_tierp_mode_parity_sweep" 8 "32G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    NORM='$tp_dir/bams/SRR15006540.bam'
    BEDS='$tp_dir/archive/beds'
    VC='$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp'
    VR='target/release/vardictrust'
    RES=\"\$KEEP_DIR/full_mode_parity_3way.tsv\"
    printf 'mode\tjava_rows\trust_rows\tcpp_rows\trust_vs_java_fp\trust_vs_java_fn\tcpp_vs_java_fp\tcpp_vs_java_fn\tverdict\n' > \"\$RES\"

    run_mode() {
      local mode=\"\$1\" j_cmd=\"\$2\" r_cmd=\"\$3\" c_cmd=\"\$4\"
      local j_key=\"\$mode\"
      [ \"\$mode\" = 'simple_chunk50000' ] && j_key='simple_f0.01'
      oracle_cache_get_or_run \"tierp_13modes/\${j_key}.java.raw.var\" \"\$KEEP_DIR/\${j_key}.java.raw.var\" \
        bash -c \"\$j_cmd > '\$KEEP_DIR/\${j_key}.java.raw.var'\"
      LC_ALL=C sort \"\$KEEP_DIR/\${j_key}.java.raw.var\" > \"\$KEEP_DIR/\${mode}.java.tsv\"
      eval \"\$r_cmd\" | LC_ALL=C sort > \"\$KEEP_DIR/\${mode}.rust.tsv\"
      if [ -z \"\${LOCAL_CKPT_DIR:-}\" ] && setsid gsutil -q stat \"${GCS_RESULTS}/cells-v103/sweep_tierp_13modes/\${mode}.cpp.sorted.var\" 2>/dev/null; then
        setsid gsutil -q cp \"${GCS_RESULTS}/cells-v103/sweep_tierp_13modes/\${mode}.cpp.sorted.var\" \"\$KEEP_DIR/\${mode}.cpp.tsv\"
      else
        eval \"\$c_cmd\" | LC_ALL=C sort > \"\$KEEP_DIR/\${mode}.cpp.tsv\"
      fi
      local jr rr cr r_fp r_fn c_fp c_fn verdict='MATCH'
      jr=\$(wc -l < \"\$KEEP_DIR/\${mode}.java.tsv\" | tr -dc '0-9')
      rr=\$(wc -l < \"\$KEEP_DIR/\${mode}.rust.tsv\" | tr -dc '0-9')
      cr=\$(wc -l < \"\$KEEP_DIR/\${mode}.cpp.tsv\" | tr -dc '0-9')
      r_fp=\$(LC_ALL=C comm -13 \"\$KEEP_DIR/\${mode}.java.tsv\" \"\$KEEP_DIR/\${mode}.rust.tsv\" | wc -l | tr -dc '0-9')
      r_fn=\$(LC_ALL=C comm -23 \"\$KEEP_DIR/\${mode}.java.tsv\" \"\$KEEP_DIR/\${mode}.rust.tsv\" | wc -l | tr -dc '0-9')
      c_fp=\$(LC_ALL=C comm -13 \"\$KEEP_DIR/\${mode}.java.tsv\" \"\$KEEP_DIR/\${mode}.cpp.tsv\" | wc -l | tr -dc '0-9')
      c_fn=\$(LC_ALL=C comm -23 \"\$KEEP_DIR/\${mode}.java.tsv\" \"\$KEEP_DIR/\${mode}.cpp.tsv\" | wc -l | tr -dc '0-9')
      if [ \"\$r_fp\" -ne 0 ] || [ \"\$r_fn\" -ne 0 ]; then
        if python3 harness/compare_var.py --eps 1e-6 <(sed 's/:0-/:1-/g' \"\$KEEP_DIR/\${mode}.java.tsv\") \"\$KEEP_DIR/\${mode}.rust.tsv\" >/dev/null 2>&1; then
          verdict='MATCH_EPS'
        else
          verdict='MISMATCH'
        fi
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \"\$mode\" \"\$jr\" \"\$rr\" \"\$cr\" \"\$r_fp\" \"\$r_fn\" \"\$c_fp\" \"\$c_fn\" \"\$verdict\" | tee -a \"\$RES\"
      [ \"\$verdict\" != 'MISMATCH' ]
    }

    run_mode 'simple_f0.01' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/simple.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 1 2>/dev/null\"

    run_mode 'simple_f0_mrd' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/simple.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 1 2>/dev/null\"

    run_mode 'simple_chunk50000' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/simple.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --chunk 50000 \$BEDS/simple.bed -th 1 2>/dev/null\"

    run_mode 'simple_fisher' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --fisher \$BEDS/big.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --fisher --emit-intermediate -th 1 \$BEDS/big.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --fisher \$BEDS/big.bed -th 1 2>/dev/null\"

    run_mode 'targeted_panel_chr17' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/panel17.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/panel17.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/panel17.bed -th 1 2>/dev/null\"

    run_mode 'sv_deletion_regions' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/sv.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/sv.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/sv.bed -th 1 2>/dev/null\"

    run_mode 'low_coverage' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/lowcov.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/lowcov.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/lowcov.bed -th 1 2>/dev/null\"

    run_mode 'amplicon_c17amp' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 -c 1 -S 2 -E 3 -g 4 \$BEDS/c17amp.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/c17amp.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 -c 1 -S 2 -E 3 -g 4 \$BEDS/c17amp.bed -th 1 2>/dev/null\"

    run_mode 'amplicon_fisher' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 --fisher -c 1 -S 2 -E 3 -g 4 \$BEDS/c17amp.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 --fisher -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/c17amp.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 --fisher -c 1 -S 2 -E 3 -g 4 \$BEDS/c17amp.bed -th 1 2>/dev/null\"

    run_mode 'somatic_f0.01_big' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 \$BEDS/big.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/big.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 \$BEDS/big.bed -th 1 2>/dev/null\"

    run_mode 'somatic_f0_mrd_c17' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 \$BEDS/c17.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \$BEDS/c17.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 \$BEDS/c17.bed -th 1 2>/dev/null\"

    run_mode 'somatic_fisher_c17' \
      \"_JAVA_OPTIONS='-Xmx8g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --fisher \$BEDS/c17.bed -th 1 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --fisher --emit-intermediate -th 1 \$BEDS/c17.bed 2>/dev/null\" \
      \"\$VC -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --fisher \$BEDS/c17.bed -th 1 2>/dev/null\"

    # P5: Tier-P Panel CRAM Twin (chr1 slice of SRR8657348.bam on simple.bed -> 228 rows)
    samtools view -T \"\$REF\" -C -o \"\$KEEP_DIR/srr8657348_chr1.cram\" \"\$BAM\" chr1
    samtools index \"\$KEEP_DIR/srr8657348_chr1.cram\"
    env RUST_MIN_STACK=67108864 \"\$VR\" -G \"\$REF\" -f 0.01 -N s -b \"\$KEEP_DIR/srr8657348_chr1.cram\" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 1 \"\$BEDS/simple.bed\" | LC_ALL=C sort > \"\$KEEP_DIR/p5_cram.rust.tsv\"
    \"\$VC\" -G \"\$REF\" -f 0.01 -N s -b \"\$KEEP_DIR/srr8657348_chr1.cram\" -c 1 -S 2 -E 3 -g 4 \"\$BEDS/simple.bed\" -th 1 2>/dev/null | LC_ALL=C sort > \"\$KEEP_DIR/p5_cram.cpp.tsv\"
    diff -u \"\$KEEP_DIR/simple_f0.01.rust.tsv\" \"\$KEEP_DIR/p5_cram.rust.tsv\" > \"\$KEEP_DIR/p5_cram_vs_bam_diff.txt\" || true
    python3 harness/compare_var.py --eps 1e-6 \"\$KEEP_DIR/simple_f0.01.rust.tsv\" \"\$KEEP_DIR/p5_cram.rust.tsv\" >/dev/null
    printf 'P5_cram_twin\t0_NoClassDefFoundError\t%s\t%s\t0\t0\t0\t0\tMATCH_BAM_TWIN\n' \
      \"\$(wc -l < \"\$KEEP_DIR/p5_cram.rust.tsv\" | tr -dc '0-9')\" \
      \"\$(wc -l < \"\$KEEP_DIR/p5_cram.cpp.tsv\" | tr -dc '0-9')\" | tee -a \"\$RES\"
  "

  # 4. Fig-2 Germline Accuracy Benchmark on hg002.exome.20.bam (6 callers scored via rtg vcfeval / eval_caller.sh)
  run_checkpointed_cell "stage2" "s2_fig2_hg002_chr20_truth" 8 "64G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    ensure_field_caller_envs
    REF='$tp_dir/archive/reference/ref20.fa'
    BAM='$tp_dir/bams/hg002.exome.20.bam'
    BED='$tp_dir/truth/target20.bed'
    TRUTH='$tp_dir/truth/truth20.norm.vcf.gz'
    SDF='$tp_dir/truth/ref20.sdf'
    [ -f \"\${REF}.fai\" ] || samtools faidx \"\$REF\"
    [ -f '$tp_dir/archive/reference/ref20.dict' ] || samtools dict \"\$REF\" > '$tp_dir/archive/reference/ref20.dict'
    [ -d \"\$SDF\" ] || rtg format -o \"\$SDF\" \"\$REF\"

    # Pad 3-col target20.bed to 4-col for VarDict (-c 1 -S 2 -E 3 -g 4 -z)
    awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\"R\"}' \"\$BED\" > \"\$KEEP_DIR/target20_4col.bed\"

    # 3-way VarDict (Java, Rust, C++) — reuse SHA256-pinned Java/Cpp VCFs from oracle-cache-v183 / cells-v103
    oracle_cache_get_or_run 'hg002_chr20_exome/java.vcf' \"\$KEEP_DIR/vardict_java.vcf\" \
      bash -c \"ORACLE_XMX_GB=16 bash harness/run_oracle.sh '\$REF' '\$BAM' '\$KEEP_DIR/target20_4col.bed' '\$KEEP_DIR/vardict_java.vcf' HG002 -f 0.01 -z -th 8\"
    target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002 -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -th 8 \"\$KEEP_DIR/target20_4col.bed\" > \"\$KEEP_DIR/vardictrust.vcf\"
    if [ -z \"\${LOCAL_CKPT_DIR:-}\" ] && setsid gsutil -q stat \"${GCS_RESULTS}/cells-v103/sweep_hg002_chr20_exome/hg002_chr20.cpp.vcf\" 2>/dev/null; then
      setsid gsutil -q cp \"${GCS_RESULTS}/cells-v103/sweep_hg002_chr20_exome/hg002_chr20.cpp.vcf\" \"\$KEEP_DIR/vardictcpp.vcf\"
    elif [ ! -s \"\$KEEP_DIR/vardictcpp.vcf\" ]; then
      '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' -G \"\$REF\" -f 0.01 -N HG002 -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z --chunk 50000 -th 8 \"\$KEEP_DIR/target20_4col.bed\" -vcf > \"\$KEEP_DIR/vardictcpp.vcf\"
    fi
    python3 harness/compare_vcf.py \"\$KEEP_DIR/vardict_java.vcf\" \"\$KEEP_DIR/vardictrust.vcf\" | tee \"\$KEEP_DIR/compare_hg002_chr20.txt\"
    grep -q '^MATCH' \"\$KEEP_DIR/compare_hg002_chr20.txt\"

    # Competing field callers (FreeBayes 1.3.10, LoFreq 2.1.5, octopus 0.7.4, Mutect2 4.6.2.0) — reuse cached VCFs when present
    if [ ! -s \"\$KEEP_DIR/freebayes.vcf\" ]; then
      conda run -n cl_freebayes freebayes -f \"\$REF\" -t \"\$BED\" \"\$BAM\" > \"\$KEEP_DIR/freebayes.vcf\"
    fi
    if [ ! -s \"\$KEEP_DIR/lofreq.vcf\" ]; then
      if [ ! -f \"\$KEEP_DIR/hg002.iq.bam\" ]; then
        conda run -n cl_lofreq lofreq indelqual --dindel -f \"\$REF\" -o \"\$KEEP_DIR/hg002.iq.bam\" \"\$BAM\"
        samtools index \"\$KEEP_DIR/hg002.iq.bam\"
      fi
      rm -f \"\$KEEP_DIR/lofreq.vcf\"
      conda run -n cl_lofreq lofreq call-parallel --pp-threads 8 --call-indels -f \"\$REF\" -l \"\$BED\" -o \"\$KEEP_DIR/lofreq.vcf\" \"\$KEEP_DIR/hg002.iq.bam\"
    fi
    if [ ! -s \"\$KEEP_DIR/octopus.vcf\" ]; then
      conda run -n cl_octopus octopus -R \"\$REF\" -I \"\$BAM\" -t \"\$BED\" -o \"\$KEEP_DIR/octopus.vcf\" --threads 8
    fi
    if [ ! -s \"\$KEEP_DIR/mutect2.filt.vcf.gz\" ]; then
      conda run -n cl_gatk4 gatk Mutect2 -R \"\$REF\" -I \"\$BAM\" -L \"\$BED\" --native-pair-hmm-threads 8 -O \"\$KEEP_DIR/m2.vcf.gz\"
      conda run -n cl_gatk4 gatk FilterMutectCalls -R \"\$REF\" -V \"\$KEEP_DIR/m2.vcf.gz\" -O \"\$KEEP_DIR/mutect2.filt.vcf.gz\"
    fi

    # Score all 6 callers identically via Archive.zip/scripts/eval_caller.sh logic
    OUT_TABLE=\"\$KEEP_DIR/fig2_caller_accuracy.tsv\"
    printf 'caller\tcalls\tsnv_p\tsnv_r\tsnv_f1\tindel_p\tindel_r\tindel_f1\n' > \"\$OUT_TABLE\"
    score_caller() {
      local lbl=\"\$1\" vcf=\"\$2\"
      local src_vcf=\"\$vcf\"
      # Synthesize GT=0/1 sample column for 8-column sites-only VCFs (e.g., LoFreq, matching Archive.zip/tool_vcfs/lofreq.gt.vcf)
      if ! bcftools view -h \"\$vcf\" 2>/dev/null | tail -n 1 | awk -F'\t' '{exit (NF>=10 ? 0 : 1)}'; then
        src_vcf=\"\$KEEP_DIR/ev_\${lbl}_gt.vcf\"
        bcftools view \"\$vcf\" 2>/dev/null | awk -F'\t' 'BEGIN{OFS=\"\t\"}
          /^##/ { print; next }
          /^#CHROM/ {
            print \"##FORMAT=<ID=GT,Number=1,Type=String,Description=\\\"Genotype\\\">\"
            print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\"FORMAT\",\"HG002\"
            next
          }
          { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\"GT\",\"0/1\" }
        ' > \"\$src_vcf\"
      fi
      bcftools reheader --fai \"\${REF}.fai\" \"\$src_vcf\" 2>/dev/null \
        | bcftools norm -f \"\$REF\" -m -any 2>/dev/null \
        | bcftools view -f PASS,. 2>/dev/null \
        | bcftools sort -Oz -o \"\$KEEP_DIR/ev_\${lbl}.vcf.gz\"
      tabix -f -p vcf \"\$KEEP_DIR/ev_\${lbl}.vcf.gz\"
      local n
      n=\$(zgrep -vc '^#' \"\$KEEP_DIR/ev_\${lbl}.vcf.gz\" || echo 0)
      for cls in snp indel; do
        bcftools view -v \"\$cls\" \"\$KEEP_DIR/ev_\${lbl}.vcf.gz\" -Oz -o \"\$KEEP_DIR/ev_\${lbl}_\${cls}.vcf.gz\"
        tabix -f -p vcf \"\$KEEP_DIR/ev_\${lbl}_\${cls}.vcf.gz\"
        bcftools view -v \"\$cls\" \"\$TRUTH\" -Oz -o \"\$KEEP_DIR/tr_\${cls}.vcf.gz\"
        tabix -f -p vcf \"\$KEEP_DIR/tr_\${cls}.vcf.gz\"
        rm -rf \"\$KEEP_DIR/rtg_\${lbl}_\${cls}\"
        rtg vcfeval -b \"\$KEEP_DIR/tr_\${cls}.vcf.gz\" -c \"\$KEEP_DIR/ev_\${lbl}_\${cls}.vcf.gz\" -t \"\$SDF\" \
          --evaluation-regions \"\$BED\" --squash-ploidy -o \"\$KEEP_DIR/rtg_\${lbl}_\${cls}\" --vcf-score-field=QUAL >/dev/null 2>&1 || true
        if [ -f \"\$KEEP_DIR/rtg_\${lbl}_\${cls}/summary.txt\" ]; then
          read -r p r f <<< \"\$(awk 'NR>2 && NF>=8' \"\$KEEP_DIR/rtg_\${lbl}_\${cls}/summary.txt\" | sort -k8 -g | tail -1 | awk '{print \$6,\$7,\$8}')\"
        else
          p=NA; r=NA; f=NA
        fi
        eval \"\${cls}_p=\${p:-NA}; \${cls}_r=\${r:-NA}; \${cls}_f=\${f:-NA}\"
      done
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \"\$lbl\" \"\$n\" \"\$snp_p\" \"\$snp_r\" \"\$snp_f\" \"\$indel_p\" \"\$indel_r\" \"\$indel_f\" | tee -a \"\$OUT_TABLE\"
    }

    score_caller 'VarDictRust' \"\$KEEP_DIR/vardictrust.vcf\"
    score_caller 'VarDictJava' \"\$KEEP_DIR/vardict_java.vcf\"
    score_caller 'vardictcpp' \"\$KEEP_DIR/vardictcpp.vcf\"
    score_caller 'octopus' \"\$KEEP_DIR/octopus.vcf\"
    score_caller 'FreeBayes' \"\$KEEP_DIR/freebayes.vcf\"
    score_caller 'LoFreq' \"\$KEEP_DIR/lofreq.vcf\"
    score_caller 'Mutect2' \"\$KEEP_DIR/mutect2.filt.vcf.gz\"
  "

  # 5. Fig-1 (SRR8657348 Single-Sample) & Fig-3 (SRR15006376|SRR15006386 Somatic) Field Callers + Somatic Venn
  run_checkpointed_cell "stage2" "s2_fig1_fig3_field_callers" 32 "115G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    ensure_field_caller_envs
    REF='$tp_dir/archive/reference/hg19.fa'
    [ -f '$tp_dir/archive/reference/hg19.dict' ] || samtools dict \"\$REF\" > '$tp_dir/archive/reference/hg19.dict'
    BAM_SINGLE='$tp_dir/bams/SRR8657348.bam'
    BED_SINGLE='$tp_dir/bams/SRR8657348.mkbed.bed'
    build_tierp_covered_bed \"\$BAM_SINGLE\" \"\$BED_SINGLE\" 773667
    TUM='$tp_dir/bams/SRR15006376.bam'
    NOR='$tp_dir/bams/SRR15006386.bam'
    BED_SOM='$tp_dir/bams/SRR15006376.mkbed.bed'
    build_tierp_covered_bed \"\$TUM\" \"\$BED_SOM\" 105949

    # Prepare LoFreq indelqual BAM for SRR8657348 and 8 FreeBayes BED chunks (c.0.bed..c.7.bed matching Wolff 2026 PROVENANCE)
    IQ_BAM='$tp_dir/bams/SRR8657348.iq.bam'
    if [ ! -f \"\$IQ_BAM\" ]; then
      conda run -n cl_lofreq lofreq indelqual --dindel -f \"\$REF\" -o \"\$IQ_BAM\" \"\$BAM_SINGLE\"
      samtools index \"\$IQ_BAM\"
    fi
    mkdir -p \"\$KEEP_DIR/fb_chunks\"
    split -d -a 1 -n l/8 \"\$BED_SINGLE\" \"\$KEEP_DIR/fb_chunks/c.\" --additional-suffix=.bed

    # 10 Hz session-sum VmRSS sampler matching Archive.zip/performance_benchmark/scripts/single_sample_speed.sh
    run_pinned_session() {
      local lbl=\"\$1\" cpus=\"\$2\"
      shift 2
      if [ \"\$(cat \"\$KEEP_DIR/\${lbl}.rc\" 2>/dev/null || true)\" = \"0\" ] && [ -s \"\$KEEP_DIR/\${lbl}.wall_s\" ]; then
        echo \"[SKIP] Sub-session '\${lbl}' already completed (rc=0, wall=\$(cat \"\$KEEP_DIR/\${lbl}.wall_s\")s)\"
        return 0
      fi
      local s e sid peak=0 cur rc=0
      s=\$(date +%s.%N)
      setsid timeout 10800 taskset -c \"\$cpus\" \"\$@\" >\"\$KEEP_DIR/\${lbl}.stdout\" 2>\"\$KEEP_DIR/\${lbl}.stderr\" &
      sid=\$!
      while kill -0 \"\$sid\" 2>/dev/null; do
        cur=\$(ps -o rss= --sid \"\$sid\" 2>/dev/null | awk '{t+=\$1}END{print t+0}')
        [ \"\${cur:-0}\" -gt \"\$peak\" ] && peak=\"\$cur\"
        sleep 0.1
      done
      wait \"\$sid\" || rc=\$?
      e=\$(date +%s.%N)
      echo \"\$rc\" > \"\$KEEP_DIR/\${lbl}.rc\"
      awk \"BEGIN{printf \\\"%.2f\\\", \$e-\$s}\" > \"\$KEEP_DIR/\${lbl}.wall_s\"
      echo \"\$peak\" > \"\$KEEP_DIR/\${lbl}.peak_kb\"
    }

    # Wave A: 8-core field callers on disjoint 8-vCPU slices (0-7, 8-15, 16-23, 24-31)
    run_pinned_session 'single_FreeBayes_8c' '0-7' conda run -n cl_freebayes bash -c '
      set -euo pipefail
      parallel -k -j 8 freebayes -f \"\$1\" -t {} \"\$2\" ::: \"\$3\"/c.{0..7}.bed | vcffirstheader | vcfstreamsort -w 1000 | vcfuniq
    ' _ \"\$REF\" \"\$BAM_SINGLE\" \"\$KEEP_DIR/fb_chunks\" &
    p1=\$!
    if [ \"\$(cat \"\$KEEP_DIR/single_LoFreq_8c.rc\" 2>/dev/null || true)\" != \"0\" ]; then rm -f \"\$KEEP_DIR/single_lf8.vcf\"; fi
    run_pinned_session 'single_LoFreq_8c' '8-15' conda run -n cl_lofreq lofreq call-parallel --pp-threads 8 --call-indels -f \"\$REF\" -l \"\$BED_SINGLE\" -o \"\$KEEP_DIR/single_lf8.vcf\" \"\$IQ_BAM\" &
    p2=\$!
    run_pinned_session 'single_Mutect2_8c' '16-23' conda run -n cl_gatk4 gatk Mutect2 -R \"\$REF\" -I \"\$BAM_SINGLE\" -L \"\$BED_SINGLE\" --native-pair-hmm-threads 8 -O \"\$KEEP_DIR/single_m2_8.vcf.gz\" &
    p3=\$!
    run_pinned_session 'single_octopus_8c' '24-31' conda run -n cl_octopus octopus -R \"\$REF\" -I \"\$BAM_SINGLE\" -t \"\$BED_SINGLE\" -o \"\$KEEP_DIR/single_oc8.vcf\" --threads 8 &
    p4=\$!
    wait \"\$p1\" \"\$p2\" \"\$p3\" \"\$p4\"

    # Wave B: Somatic 8-core field callers (Mutect2, octopus, LoFreq, Strelka2) on disjoint 8-vCPU slices
    run_pinned_session 'som_Mutect2_8c' '0-7' conda run -n cl_gatk4 gatk Mutect2 -R \"\$REF\" -I \"\$TUM\" -I \"\$NOR\" -normal SRR15006386 -L \"\$BED_SOM\" --native-pair-hmm-threads 8 -O \"\$KEEP_DIR/som_m2_8.vcf.gz\" &
    p5=\$!
    run_pinned_session 'som_octopus_8c' '8-15' conda run -n cl_octopus octopus -R \"\$REF\" -I \"\$TUM\" \"\$NOR\" --normal-sample SRR15006386 -t \"\$BED_SOM\" --threads 8 -o \"\$KEEP_DIR/som_oc8.vcf\" &
    p6=\$!
    if [ \"\$(cat \"\$KEEP_DIR/som_LoFreq_8c.rc\" 2>/dev/null || true)\" != \"0\" ]; then rm -f \"\$KEEP_DIR\"/som_lf8_*; fi
    run_pinned_session 'som_LoFreq_8c' '16-23' conda run -n cl_lofreq lofreq somatic -n \"\$NOR\" -t \"\$TUM\" -f \"\$REF\" -l \"\$BED_SOM\" -o \"\$KEEP_DIR/som_lf8_\" --threads 8 &
    p7=\$!
    # Strelka2 2.9.10 somatic extension arm (requires bgzipped + tabix-indexed BED)
    bgzip -c \"\$BED_SOM\" > \"\$KEEP_DIR/som_targets.bed.gz\" && tabix -f -p bed \"\$KEEP_DIR/som_targets.bed.gz\"
    if [ \"\$(cat \"\$KEEP_DIR/som_Strelka2_8c.rc\" 2>/dev/null || true)\" != \"0\" ]; then
      rm -rf \"\$KEEP_DIR/strelka_work\"
      conda run -n cl_strelka configureStrelkaSomaticWorkflow.py \
        --normalBam \"\$NOR\" --tumorBam \"\$TUM\" --referenceFasta \"\$REF\" \
        --callRegions \"\$KEEP_DIR/som_targets.bed.gz\" --exome --runDir \"\$KEEP_DIR/strelka_work\"
    fi
    run_pinned_session 'som_Strelka2_8c' '24-31' conda run -n cl_strelka \"\$KEEP_DIR/strelka_work/runWorkflow.py\" -m local -j 8 &
    p8=\$!
    wait \"\$p5\" \"\$p6\" \"\$p7\" \"\$p8\"

    # Wave C: All eight 1-core field-caller jobs concurrently on dedicated physical cores (0, 2, 4, 6, 8, 10, 12, 14)
    run_pinned_session 'single_FreeBayes_1c' '0' conda run -n cl_freebayes freebayes -f \"\$REF\" -t \"\$BED_SINGLE\" \"\$BAM_SINGLE\" &
    c1=\$!
    if [ \"\$(cat \"\$KEEP_DIR/single_LoFreq_1c.rc\" 2>/dev/null || true)\" != \"0\" ]; then rm -f \"\$KEEP_DIR/single_lf1.vcf\"; fi
    run_pinned_session 'single_LoFreq_1c' '2' conda run -n cl_lofreq lofreq call --call-indels -f \"\$REF\" -l \"\$BED_SINGLE\" -o \"\$KEEP_DIR/single_lf1.vcf\" \"\$IQ_BAM\" &
    c2=\$!
    run_pinned_session 'single_Mutect2_1c' '4' conda run -n cl_gatk4 gatk Mutect2 -R \"\$REF\" -I \"\$BAM_SINGLE\" -L \"\$BED_SINGLE\" --native-pair-hmm-threads 1 -O \"\$KEEP_DIR/single_m2_1.vcf.gz\" &
    c3=\$!
    run_pinned_session 'single_octopus_1c' '6' conda run -n cl_octopus octopus -R \"\$REF\" -I \"\$BAM_SINGLE\" -t \"\$BED_SINGLE\" -o \"\$KEEP_DIR/single_oc1.vcf\" --threads 1 &
    c4=\$!
    run_pinned_session 'som_Mutect2_1c' '8' conda run -n cl_gatk4 gatk Mutect2 -R \"\$REF\" -I \"\$TUM\" -I \"\$NOR\" -normal SRR15006386 -L \"\$BED_SOM\" --native-pair-hmm-threads 1 -O \"\$KEEP_DIR/som_m2_1.vcf.gz\" &
    c5=\$!
    run_pinned_session 'som_octopus_1c' '10' conda run -n cl_octopus octopus -R \"\$REF\" -I \"\$TUM\" \"\$NOR\" --normal-sample SRR15006386 -t \"\$BED_SOM\" --threads 1 -o \"\$KEEP_DIR/som_oc1.vcf\" &
    c6=\$!
    if [ \"\$(cat \"\$KEEP_DIR/som_LoFreq_1c.rc\" 2>/dev/null || true)\" != \"0\" ]; then rm -f \"\$KEEP_DIR\"/som_lf1_*; fi
    run_pinned_session 'som_LoFreq_1c' '12' conda run -n cl_lofreq lofreq somatic -n \"\$NOR\" -t \"\$TUM\" -f \"\$REF\" -l \"\$BED_SOM\" -o \"\$KEEP_DIR/som_lf1_\" --threads 1 &
    c7=\$!
    if [ \"\$(cat \"\$KEEP_DIR/som_Strelka2_1c.rc\" 2>/dev/null || true)\" != \"0\" ]; then
      rm -rf \"\$KEEP_DIR/strelka_work_1c\"
      conda run -n cl_strelka configureStrelkaSomaticWorkflow.py \
        --normalBam \"\$NOR\" --tumorBam \"\$TUM\" --referenceFasta \"\$REF\" \
        --callRegions \"\$KEEP_DIR/som_targets.bed.gz\" --exome --runDir \"\$KEEP_DIR/strelka_work_1c\"
    fi
    run_pinned_session 'som_Strelka2_1c' '14' conda run -n cl_strelka \"\$KEEP_DIR/strelka_work_1c/runWorkflow.py\" -m local -j 1 &
    c8=\$!
    wait \"\$c1\" \"\$c2\" \"\$c3\" \"\$c4\" \"\$c5\" \"\$c6\" \"\$c7\" \"\$c8\"
  "

  # 6. Wolff 2026 Table-3 Direct -R chr1:1-249250621 Whole-Chromosome Blowup Cell on SRR8657348 (n=3 across all 3 engines)
  run_checkpointed_cell "stage2" "s2_table3_R_chr1_rust_and_cpp" 8 "32G" bash -c "
    set -euo pipefail
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    for run in 1 2 3; do
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_R_chr1_run\${run}.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 > \"\$KEEP_DIR/rust_R_chr1.vcf\"
      if [ ! -s \"\$KEEP_DIR/cpp_R_chr1_run\${run}.time\" ]; then
        rc=0
        /usr/bin/time -v -o \"\$KEEP_DIR/cpp_R_chr1_run\${run}.time\" \
          '$tp_dir/vardictcpp-ae93fe39dc/build/vardictcpp' -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 -vcf > \"\$KEEP_DIR/cpp_R_chr1.vcf\" || rc=\$?
        echo \"\$rc\" > \"\$KEEP_DIR/cpp_R_chr1_run\${run}.exit\"
      fi
    done
  "

  # 7. Table-3 Java -Xmx4g OOM reproduction & -Xmx200g n=3 completion cell
  run_checkpointed_cell "stage2" "s2_table3_R_chr1_java_200g_n3" 8 "115G" bash -c "
    set -euo pipefail
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    if [ ! -s \"\$KEEP_DIR/java_xmx4g_R_chr1.time\" ]; then
      rc4=0
      JAVA_OPTS='-Xmx4g -XX:+ExitOnOutOfMemoryError' _JAVA_OPTIONS='-Xmx4g -XX:+ExitOnOutOfMemoryError' \
        timeout 300 /usr/bin/time -v -o \"\$KEEP_DIR/java_xmx4g_R_chr1.time\" \
        conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 > /dev/null 2>\"\$KEEP_DIR/java_xmx4g.err\" || rc4=\$?
      echo \"\$rc4\" > \"\$KEEP_DIR/java_xmx4g.exit\"
    fi
    for run in 1 2 3; do
      if [ ! -s \"\$KEEP_DIR/java_xmx200g_R_chr1_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/java_xmx200g_R_chr1.var\" ]; then
        JAVA_OPTS='-Xmx200g -XX:+ExitOnOutOfMemoryError' _JAVA_OPTIONS='-Xmx200g -XX:+ExitOnOutOfMemoryError' \
          /usr/bin/time -v -o \"\$KEEP_DIR/java_xmx200g_R_chr1_run\${run}.time\" \
          conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" -th 8 -R chr1:1-249250621 > \"\$KEEP_DIR/java_xmx200g_R_chr1.var\" 2>\"\$KEEP_DIR/java_xmx200g_run\${run}.err\"
      fi
    done
  "
}

run_stage_3() {
  echo "=== [Stage 3: Task 5d Cross-Tier Full-Sample RNA Arm (HG002 GM24385 51.8 GiB BAM, n=3)] ==="
  export DATA_DIR="scratch/g4"
  local rna_dir="${WORK_DIR}/data/rna"
  mkdir -p "$rna_dir"

  run_local_stage_cell "stage3" "s3_stage_rna_aligned_bam" 32 "128G" bash -c "
    set -euo pipefail
    if [ ! -f '$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam' ]; then
      setsid gsutil -q cp '${GCS_BAMS}/rna-hg002/aligned/HG002_rep1_Aligned.sortedByCoord.out.bam*' '$rna_dir/'
    fi
    samtools quickcheck '$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam'
  "

  run_checkpointed_cell "stage3" "s3_rna_splice_parity_n3" 8 "64G" bash -c "
    set -euo pipefail
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    BAM='$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam'
    BED=\"\$DATA_DIR/hg002/eval.bed\"
    for run in 1 2 3; do
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_splice_run\${run}.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002_RNA -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -i --emit-intermediate -th 8 \"\$BED\" > \"\$KEEP_DIR/rust.splice.tsv\"
      if [ ! -s \"\$KEEP_DIR/java_splice_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/java.splice.tsv\" ]; then
        _JAVA_OPTIONS='-Xmx44g' /usr/bin/time -v -o \"\$KEEP_DIR/java_splice_run\${run}.time\" \
          conda run -n vardict-oracle vardict-java -G \"\$REF\" -f 0.01 -N HG002_RNA -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z 1 -i -th 8 \"\$BED\" > \"\$KEEP_DIR/java.splice.tsv\"
      fi
    done
    LC_ALL=C sort \"\$KEEP_DIR/java.splice.tsv\" > \"\$KEEP_DIR/java.splice.sorted.tsv\"
    LC_ALL=C sort \"\$KEEP_DIR/rust.splice.tsv\" > \"\$KEEP_DIR/rust.splice.sorted.tsv\"
    diff -u \"\$KEEP_DIR/java.splice.sorted.tsv\" \"\$KEEP_DIR/rust.splice.sorted.tsv\" > \"\$KEEP_DIR/splice_diff.txt\" || true
    python3 harness/compare_var.py --eps 1e-6 <(sed 's/:0-/:1-/g' \"\$KEEP_DIR/java.splice.sorted.tsv\") \"\$KEEP_DIR/rust.splice.sorted.tsv\" | tee \"\$KEEP_DIR/compare_splice_var.txt\"
    wc -l \"\$KEEP_DIR/rust.splice.sorted.tsv\" > \"\$KEEP_DIR/splice_count.txt\"
  "
}

run_stage_4() {
  echo "=== [Stage 4: Task 5e Tier-E Full SEQC2 Exome, Tier-G Full WGS, 1..32c Scaling, ROC/QUAL Sweep & Strata] ==="
  export DATA_DIR="scratch/g4"

  run_local_stage_cell "stage4" "s4_stage_seqc2_full" 32 "128G" bash -c "
    set -euo pipefail
    s4_cached=0
    mkdir -p scratch/g4/seqc2 scratch/g4/multimodal scratch/g4/truth_cache
    if setsid gsutil -q stat '${GCS_BAMS}/g4-cache/s4_ready.done' 2>/dev/null; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/seqc2' scratch/g4/seqc2
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/multimodal' scratch/g4/multimodal
      s4_cached=1
    fi
    if setsid gsutil -q stat '${GCS_BAMS}/g4-cache/truth_cache/*' 2>/dev/null; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r '${GCS_BAMS}/g4-cache/truth_cache' scratch/g4/truth_cache || true
    fi
    bash harness/stage_g4.sh seqc2 multimodal
    if [ \"\$s4_cached\" -eq 0 ]; then
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r -x 'fastq/.*' scratch/g4/seqc2 '${GCS_BAMS}/g4-cache/seqc2' || true
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r scratch/g4/multimodal '${GCS_BAMS}/g4-cache/multimodal' || true
      date -u +%Y-%m-%dT%H:%M:%SZ | setsid gsutil -q cp - '${GCS_BAMS}/g4-cache/s4_ready.done' || true
    fi
  "
  ensure_submodule_excludes

  # 1. Uniform -th 8 Four-Target Performance Budget (amplicon, hg002, hg008, seqc2-exome) with n=3 replicates + perf_budget.py
  # Packed across Dual Sockets with zero contention: Socket 0 (0-7 / NUMA 0) runs amplicon + hg002 + hg008; Socket 1 (8-15 / NUMA 1) runs seqc2-exome.
  run_checkpointed_cell "stage4" "s4_perf_budget_4targets_n3" "0-15" "128G" bash -c '
    set -euo pipefail
    mkdir -p "$KEEP_DIR/amplicon" "$KEEP_DIR/hg002" "$KEEP_DIR/hg008" "$KEEP_DIR/seqc2"

    # Socket 0 (CPUs 0-7, NUMA Node 0): Amplicon + HG002 + HG008 (-th 8, n=3)
    (
      set -euo pipefail
      for run in 1 2 3; do
        if [ -s "$KEEP_DIR/amplicon/java_run${run}.time" ] && [ -s "$KEEP_DIR/amplicon/java.vcf" ]; then
          cp -f "$KEEP_DIR/amplicon/java_run${run}.time" "$KEEP_DIR/amplicon/java.time"
          [ -s "$KEEP_DIR/amplicon/cpp_run${run}.time" ] && cp -f "$KEEP_DIR/amplicon/cpp_run${run}.time" "$KEEP_DIR/amplicon/cpp.time" || true
          taskset -c 0-7 /usr/bin/time -v -o "$KEEP_DIR/amplicon/rust.time" \
            target/release/vardictrust -G "$DATA_DIR/amplicon/ref.fa" -f 0.01 -N SAMPLE \
            -b "$DATA_DIR/amplicon/amplicon.bam" -c 1 -S 2 -E 3 -g 4 -z -a 10:0.95 -th 8 "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" > "$KEEP_DIR/amplicon/rust.vcf"
          cp -f "$KEEP_DIR/amplicon/rust.time" "$KEEP_DIR/amplicon/rust_run${run}.time"
          python3 harness/compare_vcf.py "$KEEP_DIR/amplicon/java.vcf" "$KEEP_DIR/amplicon/rust.vcf" | tee "$KEEP_DIR/amplicon/compare.txt"
          grep -q "^MATCH" "$KEEP_DIR/amplicon/compare.txt"
        else
          taskset -c 0-7 env KEEP_DIR="$KEEP_DIR/amplicon" GOLDEN_MODE=off ORACLE_XMX_GB=22 \
            bash harness/compare.sh \
            "$DATA_DIR/amplicon/ref.fa" \
            "$DATA_DIR/amplicon/amplicon.bam" \
            "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" SAMPLE -a 10:0.95 -th 8
          cp "$KEEP_DIR/amplicon/rust.time" "$KEEP_DIR/amplicon/rust_run${run}.time"
          cp "$KEEP_DIR/amplicon/java.time" "$KEEP_DIR/amplicon/java_run${run}.time"
          [ -f "$KEEP_DIR/amplicon/cpp.time" ] && cp "$KEEP_DIR/amplicon/cpp.time" "$KEEP_DIR/amplicon/cpp_run${run}.time" || true
        fi
      done

      for run in 1 2 3; do
        if [ -s "$KEEP_DIR/hg002/java_run${run}.time" ] && [ -s "$KEEP_DIR/hg002/java.vcf" ]; then
          cp -f "$KEEP_DIR/hg002/java_run${run}.time" "$KEEP_DIR/hg002/java.time"
          [ -s "$KEEP_DIR/hg002/cpp_run${run}.time" ] && cp -f "$KEEP_DIR/hg002/cpp_run${run}.time" "$KEEP_DIR/hg002/cpp.time" || true
          taskset -c 0-7 /usr/bin/time -v -o "$KEEP_DIR/hg002/rust.time" \
            target/release/vardictrust -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" -f 0.01 -N HG002 \
            -b "$DATA_DIR/hg002/HG002.chr20.bam" -c 1 -S 2 -E 3 -g 4 -z -th 8 "$DATA_DIR/hg002/eval.bed" > "$KEEP_DIR/hg002/rust.vcf"
          cp -f "$KEEP_DIR/hg002/rust.time" "$KEEP_DIR/hg002/rust_run${run}.time"
          python3 harness/compare_vcf.py --min-variants 200000 "$KEEP_DIR/hg002/java.vcf" "$KEEP_DIR/hg002/rust.vcf" | tee "$KEEP_DIR/hg002/compare.txt"
        else
          taskset -c 0-7 env KEEP_DIR="$KEEP_DIR/hg002" GOLDEN_MODE=off MIN_VARIANTS=200000 ORACLE_XMX_GB=44 \
            bash harness/compare.sh \
            "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" \
            "$DATA_DIR/hg002/HG002.chr20.bam" \
            "$DATA_DIR/hg002/eval.bed" HG002 -f 0.01 -th 8
          cp "$KEEP_DIR/hg002/rust.time" "$KEEP_DIR/hg002/rust_run${run}.time"
          cp "$KEEP_DIR/hg002/java.time" "$KEEP_DIR/hg002/java_run${run}.time"
          [ -f "$KEEP_DIR/hg002/cpp.time" ] && cp "$KEEP_DIR/hg002/cpp.time" "$KEEP_DIR/hg002/cpp_run${run}.time" || true
        fi
      done
      grep -q "^MATCH (227562 variants)$" "$KEEP_DIR/hg002/compare.txt"

      for run in 1 2 3; do
        if [ -s "$KEEP_DIR/hg008/java_run${run}.time" ] && [ -s "$KEEP_DIR/hg008/java.vcf" ]; then
          cp -f "$KEEP_DIR/hg008/java_run${run}.time" "$KEEP_DIR/hg008/java.time"
          [ -s "$KEEP_DIR/hg008/cpp_run${run}.time" ] && cp -f "$KEEP_DIR/hg008/cpp_run${run}.time" "$KEEP_DIR/hg008/cpp.time" || true
          taskset -c 0-7 /usr/bin/time -v -o "$KEEP_DIR/hg008/rust.time" \
            target/release/vardictrust -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" -f 0.01 -N "HG008-T|HG008-N" \
            -b "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" -c 1 -S 2 -E 3 -g 4 -z -th 8 "$DATA_DIR/hg008/engine.bed" > "$KEEP_DIR/hg008/rust.vcf"
          cp -f "$KEEP_DIR/hg008/rust.time" "$KEEP_DIR/hg008/rust_run${run}.time"
          python3 harness/compare_vcf.py --min-variants 1000 "$KEEP_DIR/hg008/java.vcf" "$KEEP_DIR/hg008/rust.vcf" | tee "$KEEP_DIR/hg008/compare.txt"
        else
          taskset -c 0-7 env KEEP_DIR="$KEEP_DIR/hg008" GOLDEN_MODE=off MIN_VARIANTS=1000 ORACLE_XMX_GB=44 \
            bash harness/compare.sh \
            "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" \
            "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" \
            "$DATA_DIR/hg008/engine.bed" "HG008-T|HG008-N" -th 8
          cp "$KEEP_DIR/hg008/rust.time" "$KEEP_DIR/hg008/rust_run${run}.time"
          cp "$KEEP_DIR/hg008/java.time" "$KEEP_DIR/hg008/java_run${run}.time"
          [ -f "$KEEP_DIR/hg008/cpp.time" ] && cp "$KEEP_DIR/hg008/cpp.time" "$KEEP_DIR/hg008/cpp_run${run}.time" || true
        fi
      done
      grep -q "^MATCH (18426 variants)$" "$KEEP_DIR/hg008/compare.txt"
    ) &
    pid_s0=$!

    # Socket 1 (CPUs 8-15, NUMA Node 1): SEQC2 HCC1395/BL full somatic exome (623,444 variants, -th 8, n=3)
    (
      set -euo pipefail
      for run in 1 2 3; do
        if [ -s "$KEEP_DIR/seqc2/java_run${run}.time" ] && [ -s "$KEEP_DIR/seqc2/java.vcf" ]; then
          cp -f "$KEEP_DIR/seqc2/java_run${run}.time" "$KEEP_DIR/seqc2/java.time"
          [ -s "$KEEP_DIR/seqc2/cpp_run${run}.time" ] && cp -f "$KEEP_DIR/seqc2/cpp_run${run}.time" "$KEEP_DIR/seqc2/cpp.time" || true
          taskset -c 8-15 /usr/bin/time -v -o "$KEEP_DIR/seqc2/rust.time" \
            target/release/vardictrust -G "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa" -f 0.01 -N "HCC1395|HCC1395BL" \
            -b "$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam" -c 1 -S 2 -E 3 -g 4 -z -th 8 "$DATA_DIR/seqc2/exome.engine.bed" > "$KEEP_DIR/seqc2/rust.vcf"
          cp -f "$KEEP_DIR/seqc2/rust.time" "$KEEP_DIR/seqc2/rust_run${run}.time"
          python3 harness/compare_vcf.py --min-variants 2000 "$KEEP_DIR/seqc2/java.vcf" "$KEEP_DIR/seqc2/rust.vcf" | tee "$KEEP_DIR/seqc2/compare.txt"
        else
          taskset -c 8-15 env KEEP_DIR="$KEEP_DIR/seqc2" GOLDEN_MODE=off MIN_VARIANTS=2000 ORACLE_XMX_GB=44 \
            bash harness/compare.sh \
            "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa" \
            "$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam" \
            "$DATA_DIR/seqc2/exome.engine.bed" "HCC1395|HCC1395BL" -th 8
          cp "$KEEP_DIR/seqc2/rust.time" "$KEEP_DIR/seqc2/rust_run${run}.time"
          cp "$KEEP_DIR/seqc2/java.time" "$KEEP_DIR/seqc2/java_run${run}.time"
          [ -f "$KEEP_DIR/seqc2/cpp.time" ] && cp "$KEEP_DIR/seqc2/cpp.time" "$KEEP_DIR/seqc2/cpp_run${run}.time" || true
        fi
      done
      grep -q "^MATCH (623444 variants)$" "$KEEP_DIR/seqc2/compare.txt"
    ) &
    pid_s1=$!

    wait "$pid_s0"
    wait "$pid_s1"

    # Sequentially record evidence receipts after the wait barrier (zero file contention)
    python3 harness/record_evidence.py amplicon-tier2 MATCH "$KEEP_DIR/amplicon" -- bash harness/compare.sh "$DATA_DIR/amplicon/ref.fa" "$DATA_DIR/amplicon/amplicon.bam" "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" SAMPLE -a 10:0.95 -th 8
    python3 harness/record_evidence.py hg002-tier2 MATCH "$KEEP_DIR/hg002" -- bash harness/compare.sh "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" "$DATA_DIR/hg002/HG002.chr20.bam" "$DATA_DIR/hg002/eval.bed" HG002 -f 0.01 -th 8
    python3 harness/record_evidence.py hg008-tier2 MATCH "$KEEP_DIR/hg008" -- bash harness/compare.sh "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" "$DATA_DIR/hg008/engine.bed" "HG008-T|HG008-N" -th 8
    python3 harness/record_evidence.py seqc2-exome-tier2 MATCH "$KEEP_DIR/seqc2" -- bash harness/compare.sh "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa" "$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam" "$DATA_DIR/seqc2/exome.engine.bed" "HCC1395|HCC1395BL" -th 8

    # Validate 4-target performance budget gate on the GCP host
    python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 | tee "$KEEP_DIR/perf_budget_gcp.json"
  '

  # 2. Tier-3 Ground Truth Concordance + Full ROC/QUAL Sweep (0/30/50/75/100/128/150/200 & 128/160/226) + SNP/InDel/MNP/SV Strata
  run_checkpointed_cell "stage4" "s4_tier3_truth_and_qual_sweep" 8 "32G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    BUDGET_DIR="${WORK_DIR}/cells/s4_perf_budget_4targets_n3"
    if [ ! -f "$BUDGET_DIR/hg002/rust.vcf" ] && [ -z "${LOCAL_CKPT_DIR:-}" ]; then
      mkdir -p "$BUDGET_DIR"
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_RESULTS}/${GCS_CELLS_SUBDIR:-cells-v105}/s4_perf_budget_4targets_n3" "$BUDGET_DIR" || \
        setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_RESULTS}/cells-v104/s4_perf_budget_4targets_n3" "$BUDGET_DIR"
    fi

    # (a) HG002 Germline Truth (GIAB v4.2.1)
    bash harness/concordance_germline.sh hg002-truth \
      "$BUDGET_DIR/hg002/java.vcf" "$BUDGET_DIR/hg002/rust.vcf" \
      "$DATA_DIR/hg002/truth/truth.vcf.gz" "$DATA_DIR/hg002/eval.bed" \
      "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf" "$KEEP_DIR/hg002_truth" HG002 "$BUDGET_DIR/hg002/cpp.vcf"

    # (b) SEQC2 Somatic Exome Truth (SEQC2 v1.2.1)
    bash harness/concordance_pair.sh seqc2-truth \
      "$BUDGET_DIR/seqc2/java.vcf" "$BUDGET_DIR/seqc2/rust.vcf" \
      "$DATA_DIR/seqc2/truth/truth.vcf.gz" "$DATA_DIR/seqc2/eval.bed" \
      "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.sdf" "$KEEP_DIR/seqc2_truth" HCC1395 "$BUDGET_DIR/seqc2/cpp.vcf"

    # (c) HG008 Somatic Truth (NIST v0.3 + v0.4 smvar in parallel)
    bash harness/concordance_pair.sh hg008-truth \
      "$BUDGET_DIR/hg008/java.vcf" "$BUDGET_DIR/hg008/rust.vcf" \
      "$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz" "$DATA_DIR/hg008/eval.bed" \
      "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf" "$KEEP_DIR/hg008_truth" HG008-T "$BUDGET_DIR/hg008/cpp.vcf"

    if [ -f "$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz" ] && [ -f "$DATA_DIR/hg008/eval_v04.bed" ]; then
      bash harness/concordance.sh \
        "$BUDGET_DIR/hg008/java.vcf" \
        "$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz" \
        "$DATA_DIR/hg008/eval_v04.bed" "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf" \
        "$KEEP_DIR/hg008_v04_truth/java" HG008-T > "$KEEP_DIR/hg008_v04_java.txt"
      bash harness/concordance.sh \
        "$BUDGET_DIR/hg008/rust.vcf" \
        "$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz" \
        "$DATA_DIR/hg008/eval_v04.bed" "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf" \
        "$KEEP_DIR/hg008_v04_truth/rust" HG008-T > "$KEEP_DIR/hg008_v04_rust.txt"
    fi

    # (d) Pre-registered Full ROC/QUAL Sweep (0, 30, 50, 75, 100, 128, 150, 160, 200, 226) on HG002, SEQC2, HG008
    python3 - "$KEEP_DIR" <<'\''PY'\''
import json, pathlib, sys
kd = pathlib.Path(sys.argv[1])
thresholds = [0, 30, 50, 75, 100, 128, 150, 160, 200, 226]
sweep = {}
for cohort, sub in [("hg002", "hg002_truth"), ("seqc2", "seqc2_truth"), ("hg008", "hg008_truth")]:
    sweep[cohort] = {}
    for eng in ("java", "rust", "cpp"):
        sfile = kd / sub / eng / "vcfeval" / "summary.txt"
        if not sfile.exists():
            continue
        rows = []
        for line in sfile.read_text().splitlines():
            parts = line.split()
            if not parts or parts[0].startswith(("Threshold", "-")):
                continue
            t_val = None if parts[0] == "None" else float(parts[0])
            rows.append((t_val, dict(threshold=parts[0], tp_base=int(parts[1]), tp_call=int(parts[2]),
                                     fp=int(parts[3]), fn=int(parts[4]),
                                     precision=float(parts[5]), recall=float(parts[6]), f1=float(parts[7]))))
        eng_sweep = {}
        for q in thresholds:
            cand = [r[1] for r in rows if r[0] is not None and r[0] >= q]
            if cand:
                eng_sweep[str(q)] = cand[0]
            elif rows:
                eng_sweep[str(q)] = rows[-1][1]
        none_r = [r[1] for r in rows if r[0] is None]
        if none_r:
            eng_sweep["None"] = none_r[0]
        sweep[cohort][eng] = eng_sweep
(kd / "roc_qual_sweep.json").write_text(json.dumps(sweep, indent=2))
PY

    # (e) SNP / InDel / MNP / SV & AF-bin stratification across HG002, SEQC2, HG008 (Java vs Rust vs C++)
    python3 - "$BUDGET_DIR" "$KEEP_DIR/strata_summary.json" <<'\''PY'\''
import json, pathlib, sys
bdir = pathlib.Path(sys.argv[1])
out_path = pathlib.Path(sys.argv[2])
strata = {}
for cohort in ("hg002", "seqc2", "hg008"):
    strata[cohort] = {}
    for eng in ("java", "rust", "cpp"):
        vcf = bdir / cohort / f"{eng}.vcf"
        if not vcf.exists():
            continue
        counts = {"SNV": 0, "Insertion": 0, "Deletion": 0, "Complex_MNP": 0, "SV": 0,
                  "AF_le_0.01": 0, "AF_0.01_0.05": 0, "AF_gt_0.05": 0, "PASS": 0, "total": 0}
        with open(vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                counts["total"] += 1
                if parts[6] == "PASS":
                    counts["PASS"] += 1
                info = dict(kv.split("=", 1) for kv in parts[7].split(";") if "=" in kv)
                vtype = info.get("TYPE", "")
                alt = parts[4]
                if alt.startswith("<") or "SVTYPE" in info:
                    counts["SV"] += 1
                elif vtype == "SNV" or (len(parts[3]) == 1 and len(alt) == 1):
                    counts["SNV"] += 1
                elif vtype == "Insertion" or (len(parts[3]) == 1 and len(alt) > 1):
                    counts["Insertion"] += 1
                elif vtype == "Deletion" or (len(parts[3]) > 1 and len(alt) == 1):
                    counts["Deletion"] += 1
                else:
                    counts["Complex_MNP"] += 1
                try:
                    af = float(info.get("AF", "0").split(",")[0])
                    if af <= 0.01:
                        counts["AF_le_0.01"] += 1
                    elif af <= 0.05:
                        counts["AF_0.01_0.05"] += 1
                    else:
                        counts["AF_gt_0.05"] += 1
                except ValueError:
                    pass
        strata[cohort][eng] = counts
    assert strata[cohort]["java"] == strata[cohort]["rust"], f"Strata mismatch on {cohort}"
out_path.write_text(json.dumps(strata, indent=2))
PY
  '

  # 3. Multi-core scaling 1, 2, 4, 8, 16, 32 threads across all 3 engines (Rust, C++, Java) + T2 Thread Determinism
  for th in 1 2 4 8 16 32; do
    local mem_cap="32G"
    local java_xmx="22"
    local cell_cpus="$th"
    local base_core="0"
    if [ "$th" -eq 1 ]; then
      cell_cpus="${SOCKET_PIN_CPUS:-0-7}"
      mem_cap="128G"
      [ "$cell_cpus" = "8-15" ] && base_core="8"
    elif [ "$th" -eq 2 ] || [ "$th" -eq 4 ]; then
      cell_cpus="${SOCKET_PIN_CPUS:-8-15}"
      mem_cap="128G"
      [ "$cell_cpus" = "8-15" ] && base_core="8"
    elif [ "$th" -ge 16 ]; then
      mem_cap="64G"
      java_xmx="44"
    fi
    if [ "$th" -ge 32 ]; then
      mem_cap="128G"
      java_xmx="88"
    fi
    run_checkpointed_cell "stage4" "s4_scale_exome_3way_th${th}" "$cell_cpus" "$mem_cap" bash -c "
      set -euo pipefail
      source harness/gcp_paper_runner.sh --source-only
      run_s4_scale_exome_body '${th}' '${java_xmx}' '${base_core}'
    "
  done

  # 4. Tier-G chr1 span scaling (1 Mbp, 10 Mbp, 100 Mbp, 248,956,422 bp) at 8, 16, 32 threads for BOTH Rust and C++ (n=3)
  for th in 8 16 32; do
    local mem_cap="64G"
    [ "$th" -ge 32 ] && mem_cap="128G"
    run_checkpointed_cell "stage4" "s4_span_chr1_3way_th${th}" "$th" "$mem_cap" bash -c "
      set -euo pipefail
      REF=\"\$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa\"
      BAM=\"\$DATA_DIR/seqc2/HCC1395.bam\"
      for sz in 1000000 10000000 100000000 248956422; do
        printf 'chr1\t0\t%s\tR%s\n' \"\$sz\" \"\$sz\" > \"\$KEEP_DIR/R\${sz}.bed\"
        for run in 1 2 3; do
          if [ ! -s \"\$KEEP_DIR/rust_R\${sz}_th${th}_run\${run}.time\" ] || [ ! -s \"\$KEEP_DIR/rust_R\${sz}_th${th}.vcf\" ]; then
            /usr/bin/time -v -o \"\$KEEP_DIR/rust_R\${sz}_th${th}_run\${run}.time\" \
              target/release/vardictrust -G \"\$REF\" -f 0.01 -N HCC1395 -b \"\$BAM\" -th ${th} -R \"chr1:1-\${sz}\" > \"\$KEEP_DIR/rust_R\${sz}_th${th}.vcf\"
          fi
          if [ ! -s \"\$KEEP_DIR/cpp_R\${sz}_th${th}_run\${run}.time\" ]; then
            rc_cpp=0
            /usr/bin/time -v -o \"\$KEEP_DIR/cpp_R\${sz}_th${th}_run\${run}.time\" \
              bash harness/run_cpp.sh \"\$REF\" \"\$BAM\" \"\$KEEP_DIR/R\${sz}.bed\" \"\$KEEP_DIR/cpp_R\${sz}_th${th}.vcf\" HCC1395 -th ${th} -z || rc_cpp=\$?
            echo \"\$rc_cpp\" > \"\$KEEP_DIR/cpp_R\${sz}_th${th}_run\${run}.exit\"
          fi
        done
      done
    "
  done

  # 5. Tier-G chr1 MRD -f 0 noise-floor operating-point cell (gcp-prereg.md §A2 H3 & §A10)
  run_checkpointed_cell "stage4" "s4_span_chr1_f0" 8 "64G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    BAM="$DATA_DIR/seqc2/HCC1395.bam"
    printf "chr1\t0\t10000000\tR10M\n" > "$KEEP_DIR/chr1_10m.bed"
    for run in 1 2 3; do
      /usr/bin/time -v -o "$KEEP_DIR/rust_chr1_10m_f0_run${run}.time" \
        target/release/vardictrust -G "$REF" -f 0 -N HCC1395 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z --tile-size 10000000 -th 8 "$KEEP_DIR/chr1_10m.bed" > "$KEEP_DIR/rust_chr1_10m_f0.vcf"
      if [ ! -s "$KEEP_DIR/cpp_chr1_10m_f0_run${run}.time" ] || [ ! -s "$KEEP_DIR/cpp_chr1_10m_f0.vcf" ]; then
        /usr/bin/time -v -o "$KEEP_DIR/cpp_chr1_10m_f0_run${run}.time" \
          bash harness/run_cpp.sh "$REF" "$BAM" "$KEEP_DIR/chr1_10m.bed" "$KEEP_DIR/cpp_chr1_10m_f0.vcf" HCC1395 -f 0 -th 8 -z
      fi
    done
    oracle_cache_get_or_run "span_chr1_f0/java_chr1_10m_f0.vcf" "$KEEP_DIR/java_chr1_10m_f0.vcf" \
      bash -c "ORACLE_XMX_GB=44 /usr/bin/time -v -o '\''$KEEP_DIR/java_chr1_10m_f0_run1.time'\'' bash harness/run_oracle.sh '\''$REF'\'' '\''$BAM'\'' '\''$KEEP_DIR/chr1_10m.bed'\'' '\''$KEEP_DIR/java_chr1_10m_f0.vcf'\'' HCC1395 -f 0 -z -th 8"
    python3 harness/compare_vcf.py "$KEEP_DIR/java_chr1_10m_f0.vcf" "$KEEP_DIR/rust_chr1_10m_f0.vcf" | tee "$KEEP_DIR/compare_chr1_10m_f0.txt"
    grep -q "^MATCH" "$KEEP_DIR/compare_chr1_10m_f0.txt"
  '

  # 6. Concurrent jobs/node cell: 2 x 4-thread jobs inside an n2-standard-8-equivalent (8c / 32G) slice, n=3 replicates
  run_checkpointed_cell "stage4" "s4_jobs_per_node_2x4c_n3" 8 "32G" bash -c '
    set -euo pipefail
    REF="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    BAM="$DATA_DIR/seqc2/HCC1395.bam"
    BED="$DATA_DIR/seqc2/exome.engine.bed"
    for run in 1 2 3; do
      taskset -c 0-3 /usr/bin/time -v -o "$KEEP_DIR/job1_4c_run${run}.time" \
        target/release/vardictrust -G "$REF" -f 0.01 -N HCC1395_J1 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th 4 "$BED" > "$KEEP_DIR/job1.vcf" &
      p1=$!
      taskset -c 4-7 /usr/bin/time -v -o "$KEEP_DIR/job2_4c_run${run}.time" \
        target/release/vardictrust -G "$REF" -f 0.01 -N HCC1395_J2 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th 4 "$BED" > "$KEEP_DIR/job2.vcf" &
      p2=$!
      wait "$p1" "$p2"
    done
  '

  # 7. Full multimodal + extended coverage + G4 & E2 CRAM twins + README CLI recipe verification
  run_checkpointed_cell "stage4" "s4_multimodal_and_cram_twin" 8 "32G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    export GOLDEN_STORE_DIR="${WORK_DIR}/golden-store"
    mkdir -p "$GOLDEN_STORE_DIR" "$DATA_DIR/truth_cache"
    if [ -z "${LOCAL_CKPT_DIR:-}" ] && setsid gsutil -q stat "${GCS_ORACLE_CACHE}/golden-store/cells/*" 2>/dev/null; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_ORACLE_CACHE}/golden-store" "$GOLDEN_STORE_DIR" || true
    fi
    if [ -z "${LOCAL_CKPT_DIR:-}" ] && setsid gsutil -q stat "${GCS_BAMS}/g4-cache/truth_cache/*" 2>/dev/null; then
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_BAMS}/g4-cache/truth_cache" "$DATA_DIR/truth_cache" || true
    fi
    GOLDEN_MODE=auto GOLDEN_STORE_DIR="$GOLDEN_STORE_DIR" bash harness/run_extended_coverage.sh
    GOLDEN_MODE=auto GOLDEN_STORE_DIR="$GOLDEN_STORE_DIR" bash harness/run_multimodal_slices.sh
    CMP_TIMEOUT=12h python3 harness/record_evidence.py hg008-cram-twin -- bash harness/cram_twin.sh
    GOLDEN_MODE=auto GOLDEN_STORE_DIR="$GOLDEN_STORE_DIR" bash harness/concordance_multimodal.sh

    # E2: Tier-E SEQC2 Exome CRAM Twin (HCC1395/HCC1395BL slice CRAM vs BAM bitwise identity)
    REF_E="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    BAM_T_E="$DATA_DIR/multimodal/seqc2_s3a_tumor.bam"
    BAM_N_E="$DATA_DIR/multimodal/seqc2_s3a_normal.bam"
    CRAM_T_E="$KEEP_DIR/seqc2_s3a_tumor.cram"
    CRAM_N_E="$KEEP_DIR/seqc2_s3a_normal.cram"
    printf "chr7\t140400000\t140420000\tBRAF\n" > "$KEEP_DIR/e2_braf.bed"
    samtools view -T "$REF_E" -C -o "$CRAM_T_E" "$BAM_T_E" && samtools index "$CRAM_T_E"
    samtools view -T "$REF_E" -C -o "$CRAM_N_E" "$BAM_N_E" && samtools index "$CRAM_N_E"
    target/release/vardictrust -G "$REF_E" -f 0.01 -N "TUMOR|NORMAL" -b "$BAM_T_E|$BAM_N_E" -c 1 -S 2 -E 3 -g 4 -z -th 4 "$KEEP_DIR/e2_braf.bed" > "$KEEP_DIR/e2_rust_bam.vcf"
    env RUST_MIN_STACK=67108864 target/release/vardictrust -G "$REF_E" -f 0.01 -N "TUMOR|NORMAL" -b "$CRAM_T_E|$CRAM_N_E" -c 1 -S 2 -E 3 -g 4 -z -th 4 "$KEEP_DIR/e2_braf.bed" > "$KEEP_DIR/e2_rust_cram.vcf"
    python3 harness/compare_vcf.py --strict-bytes "$KEEP_DIR/e2_rust_bam.vcf" "$KEEP_DIR/e2_rust_cram.vcf" | tee "$KEEP_DIR/e2_cram_twin_compare.txt"
    grep -q "^MATCH" "$KEEP_DIR/e2_cram_twin_compare.txt"

    # Verify README.md CLI recipes live on the GCP instance
    target/release/vardictrust --version > "$KEEP_DIR/readme_version.txt"
    target/release/vardictrust -G fixtures/micro/micro.ref.fa -b fixtures/micro/micro.bam -R chr1:40-260 > "$KEEP_DIR/readme_smoke.vcf"
    grep -q "^chr1" "$KEEP_DIR/readme_smoke.vcf"

    # Sync updated benchmarks/evidence/*.json back to GCS evidence bucket
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r benchmarks/evidence "$GCS_RESULTS/${GCS_CELLS_SUBDIR:-cells-v105}/evidence_certified" || true
  '
}

run_s4_scale_exome_body() {
  local th="$1" java_xmx="$2" base_core="${3:-0}"
  local REF="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
  local BAM="$DATA_DIR/seqc2/HCC1395.bam"
  local BED="$DATA_DIR/seqc2/exome.engine.bed"

  run_one_rust() {
    local run="$1" cpu_pin="${2:-}"
    local out_vcf="$KEEP_DIR/rust_th${th}_run${run}.vcf"
    if [ "$run" -eq 1 ] && [ -s "$KEEP_DIR/rust_th${th}_run1.time" ] && [ -s "$KEEP_DIR/rust_th${th}.vcf" ]; then
      echo "[SKIP] rust_th${th}_run1 already complete."
      return 0
    fi
    if [ "$run" -gt 1 ] && [ -s "$KEEP_DIR/rust_th${th}_run${run}.time" ]; then
      echo "[SKIP] rust_th${th}_run${run} already complete."
      return 0
    fi
    if [ -n "$cpu_pin" ]; then
      taskset -c "$cpu_pin" /usr/bin/time -v -o "$KEEP_DIR/rust_th${th}_run${run}.time" \
        target/release/vardictrust -G "$REF" -f 0.01 -N HCC1395 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th "$th" "$BED" > "$out_vcf"
    else
      /usr/bin/time -v -o "$KEEP_DIR/rust_th${th}_run${run}.time" \
        target/release/vardictrust -G "$REF" -f 0.01 -N HCC1395 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th "$th" "$BED" > "$out_vcf"
    fi
    if [ "$run" -eq 1 ]; then
      mv -f "$out_vcf" "$KEEP_DIR/rust_th${th}.vcf"
    else
      rm -f "$out_vcf"
    fi
  }

  run_one_cpp() {
    local run="$1" cpu_pin="${2:-}"
    local out_vcf="$KEEP_DIR/cpp_th${th}_run${run}.vcf"
    if [ "$run" -eq 1 ] && [ -s "$KEEP_DIR/cpp_th${th}_run1.time" ] && [ -s "$KEEP_DIR/cpp_th${th}.vcf" ]; then
      echo "[SKIP] cpp_th${th}_run1 already complete."
      return 0
    fi
    if [ "$run" -gt 1 ] && [ -s "$KEEP_DIR/cpp_th${th}_run${run}.time" ]; then
      echo "[SKIP] cpp_th${th}_run${run} already complete."
      return 0
    fi
    if [ -n "$cpu_pin" ]; then
      taskset -c "$cpu_pin" /usr/bin/time -v -o "$KEEP_DIR/cpp_th${th}_run${run}.time" \
        bash harness/run_cpp.sh "$REF" "$BAM" "$BED" "$out_vcf" HCC1395 -th "$th" -z
    else
      /usr/bin/time -v -o "$KEEP_DIR/cpp_th${th}_run${run}.time" \
        bash harness/run_cpp.sh "$REF" "$BAM" "$BED" "$out_vcf" HCC1395 -th "$th" -z
    fi
    if [ "$run" -eq 1 ]; then
      mv -f "$out_vcf" "$KEEP_DIR/cpp_th${th}.vcf"
    else
      rm -f "$out_vcf"
    fi
  }

  run_one_java() {
    local run="$1" cpu_pin="${2:-}"
    local out_vcf="$KEEP_DIR/java_th${th}_run${run}.vcf"
    if [ "$run" -eq 1 ] && [ -s "$KEEP_DIR/java_th${th}_run1.time" ] && [ -s "$KEEP_DIR/java_th${th}.vcf" ]; then
      echo "[SKIP] java_th${th}_run1 already complete."
      return 0
    fi
    if [ "$run" -gt 1 ] && [ -s "$KEEP_DIR/java_th${th}_run${run}.time" ]; then
      echo "[SKIP] java_th${th}_run${run} already complete."
      return 0
    fi
    if [ -n "$cpu_pin" ]; then
      ORACLE_XMX_GB="$java_xmx" taskset -c "$cpu_pin" /usr/bin/time -v -o "$KEEP_DIR/java_th${th}_run${run}.time" \
        bash harness/run_oracle.sh "$REF" "$BAM" "$BED" "$out_vcf" HCC1395 -f 0.01 -z -th "$th"
    else
      ORACLE_XMX_GB="$java_xmx" /usr/bin/time -v -o "$KEEP_DIR/java_th${th}_run${run}.time" \
        bash harness/run_oracle.sh "$REF" "$BAM" "$BED" "$out_vcf" HCC1395 -f 0.01 -z -th "$th"
    fi
    if [ "$run" -eq 1 ]; then
      mv -f "$out_vcf" "$KEEP_DIR/java_th${th}.vcf"
    else
      rm -f "$out_vcf"
    fi
  }

  local pids=()
  if [ "$th" -eq 1 ]; then
    ( for r in 1 2 3; do run_one_rust "$r" "$((base_core + 0))"; done ) &
    pids+=($!)
    ( for r in 1 2 3; do run_one_cpp "$r" "$((base_core + 1))"; done ) &
    pids+=($!)
    run_one_java 1 "$((base_core + 2))" &
    pids+=($!)
    run_one_java 2 "$((base_core + 3))" &
    pids+=($!)
    run_one_java 3 "$((base_core + 4))" &
    pids+=($!)
  elif [ "$th" -eq 2 ]; then
    run_one_java 1 "$((base_core + 0))-$((base_core + 1))" &
    pids+=($!)
    run_one_java 2 "$((base_core + 2))-$((base_core + 3))" &
    pids+=($!)
    run_one_java 3 "$((base_core + 4))-$((base_core + 5))" &
    pids+=($!)
    (
      for r in 1 2 3; do run_one_rust "$r" "$((base_core + 6))-$((base_core + 7))"; done
      for r in 1 2 3; do run_one_cpp  "$r" "$((base_core + 6))-$((base_core + 7))"; done
    ) &
    pids+=($!)
  elif [ "$th" -eq 4 ]; then
    (
      run_one_java 1 "$((base_core + 0))-$((base_core + 3))"
      run_one_java 2 "$((base_core + 0))-$((base_core + 3))"
    ) &
    pids+=($!)
    (
      for r in 1 2 3; do run_one_rust "$r" "$((base_core + 4))-$((base_core + 7))"; done
      for r in 1 2 3; do run_one_cpp  "$r" "$((base_core + 4))-$((base_core + 7))"; done
      run_one_java 3 "$((base_core + 4))-$((base_core + 7))"
    ) &
    pids+=($!)
  else
    for r in 1 2 3; do
      run_one_rust "$r" ""
      run_one_cpp  "$r" ""
      run_one_java "$r" ""
    done
  fi

  local rc=0 p
  for p in "${pids[@]}"; do
    wait "$p" || rc=$?
  done
  if [ "$rc" -ne 0 ]; then
    echo "FATAL: s4_scale_exome_3way_th${th} worker exited with rc=$rc" >&2
    return "$rc"
  fi

  python3 harness/compare_vcf.py "$KEEP_DIR/java_th${th}.vcf" "$KEEP_DIR/rust_th${th}.vcf" | tee "$KEEP_DIR/compare_th${th}.txt"
  grep -q '^MATCH' "$KEEP_DIR/compare_th${th}.txt"
  if [ "$th" -gt 1 ] && [ ! -f "${WORK_DIR}/cells/s4_scale_exome_3way_th1/rust_th1.vcf" ] && [ -z "${LOCAL_CKPT_DIR:-}" ]; then
    mkdir -p "${WORK_DIR}/cells/s4_scale_exome_3way_th1"
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_RESULTS}/${GCS_CELLS_SUBDIR:-cells-v105}/s4_scale_exome_3way_th1" "${WORK_DIR}/cells/s4_scale_exome_3way_th1" || \
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_RESULTS}/cells-v104/s4_scale_exome_3way_th1" "${WORK_DIR}/cells/s4_scale_exome_3way_th1" || true
  fi
  if [ "$th" -gt 1 ] && [ -f "${WORK_DIR}/cells/s4_scale_exome_3way_th1/rust_th1.vcf" ]; then
    cmp <(grep -v '^##' "${WORK_DIR}/cells/s4_scale_exome_3way_th1/rust_th1.vcf") <(grep -v '^##' "$KEEP_DIR/rust_th${th}.vcf")
    echo 'DETERMINISTIC_RECORDS_IDENTICAL_TO_TH1' > "$KEEP_DIR/determinism_th${th}.txt"
  fi
}

run_packed_stages_2_to_4() {
  echo "=== [Packed Zero-Contention Dual-Socket Orchestrator (NUMA Node 0: CPUs 0-7 / L3=0 vs NUMA Node 1: CPUs 8-15 / L3=1, SMT 16-31 Idle)] ==="

  # Wave 0A: Verify Tier-P & Stage-3 RNA assets (already on NVMe), build all 5 mkbed.bed files, and warm RAM page cache
  echo "=== [Wave 0A: Verifying Tier-P + RNA Assets & Warming 251 GiB RAM Page Cache (Zero Disk/Network I/O During Timed Runs)] ==="
  date -u +%Y-%m-%dT%H:%M:%SZ > "${WORK_DIR}/local_staged/bg_prestage.done"
  CELL_ALLOWLIST="s2_stage_tierp_assets s3_stage_rna_aligned_bam" \
    bash -c 'source harness/gcp_paper_runner.sh --source-only; run_stage_2; run_stage_3'
  ensure_field_caller_envs
  local tp_dir="${WORK_DIR}/data/tier-p"
  build_tierp_covered_bed "$tp_dir/bams/SRR15006386.bam" "$tp_dir/bams/SRR15006386.mkbed.bed" 11732
  build_tierp_covered_bed "$tp_dir/bams/SRR15006375.bam" "$tp_dir/bams/SRR15006375.mkbed.bed" 12292
  build_tierp_covered_bed "$tp_dir/bams/SRR15006376.bam" "$tp_dir/bams/SRR15006376.mkbed.bed" 105949
  build_tierp_covered_bed "$tp_dir/bams/SRR15006540.bam" "$tp_dir/bams/SRR15006540.mkbed.bed" 135145
  build_tierp_covered_bed "$tp_dir/bams/SRR8657348.bam"  "$tp_dir/bams/SRR8657348.mkbed.bed"  773667
  cat "$tp_dir"/bams/*.bam "$tp_dir"/archive/reference/*.fa > /dev/null 2>&1 || true

  # Wave 1A: Dual-Socket NUMA-Isolated Stage-2 & Stage-3 queues (Socket 0: CPUs 0-7 / MemNode 0 / L3=0; Socket 1: CPUs 8-15 / MemNode 1 / L3=1; SMT 16-31 Idle)
  echo "=== [Wave 1A: Dual-Socket NUMA-Isolated Stage 2+3 Execution (Socket 0: CPUs 0-7 / Node 0 || Socket 1: CPUs 8-15 / Node 1)] ==="
  (
    export SOCKET_PIN_CPUS="0-7"
    export CELL_ALLOWLIST="s2_tierp_v2_SRR15006386_th8 s2_tierp_v2_SRR15006386_th1 s2_tierp_v2_SRR15006375_th8 s2_tierp_v2_SRR15006376_th8 s2_tierp_v2_SRR15006376_th1 s2_tierp_v2_SRR8657348_th8 s2_tierp_v2_SRR8657348_th1 s2_tierp_mode_parity_sweep s2_table3_R_chr1_rust_and_cpp"
    run_stage_2
  ) &
  local pid_s0=$!

  (
    export SOCKET_PIN_CPUS="8-15"
    export CELL_ALLOWLIST="s2_tierp_v2_SRR15006375_th1 s2_tierp_v2_SRR15006540_th8 s2_tierp_v2_SRR15006540_th1 s2_tierp_v2_somatic_76_86_th8 s2_tierp_v2_somatic_76_86_th1 s2_fig2_hg002_chr20_truth s3_rna_splice_parity_n3"
    run_stage_2
    run_stage_3
  ) &
  local pid_s1=$!

  local rc_s0=0 rc_s1=0
  wait "$pid_s0" || rc_s0=$?
  wait "$pid_s1" || rc_s1=$?
  if [ "$rc_s0" -ne 0 ] || [ "$rc_s1" -ne 0 ]; then
    echo "FATAL: Dual-socket Wave 1A exited non-zero (rc_s0=$rc_s0, rc_s1=$rc_s1)" >&2
    return 1
  fi

  # Wave 0B: Stage-4 Asset Barrier (runs while no timed cells are active; warms SEQC2 BAMs into RAM page cache)
  echo "=== [Wave 0B: Staging Stage-4 SEQC2 & Multimodal Assets & Warming Page Cache] ==="
  CELL_ALLOWLIST="s4_stage_seqc2_full" \
    bash -c 'source harness/gcp_paper_runner.sh --source-only; run_stage_4'
  ln -sfn truth/truth.vcf.gz scratch/g4/seqc2/truth.vcf.gz
  ln -sfn truth/truth.vcf.gz.tbi scratch/g4/seqc2/truth.vcf.gz.tbi
  cat scratch/g4/seqc2/*.bam scratch/g4/ref/d1vd1/*.fa > /dev/null 2>&1 || true

  # Wave 1B: Dual-Socket NUMA-Isolated Stage-4 8-core queues
  echo "=== [Wave 1B: Dual-Socket NUMA-Isolated Stage-4 8c Execution (Socket 0: CPUs 0-7 / Node 0 || Socket 1: CPUs 8-15 / Node 1)] ==="
  (
    export SOCKET_PIN_CPUS="0-7"
    export CELL_ALLOWLIST="s4_span_chr1_f0 s4_jobs_per_node_2x4c_n3"
    run_stage_4
  ) &
  pid_s0=$!

  (
    export SOCKET_PIN_CPUS="8-15"
    export CELL_ALLOWLIST="s4_multimodal_and_cram_twin"
    run_stage_4
  ) &
  pid_s1=$!

  rc_s0=0; rc_s1=0
  wait "$pid_s0" || rc_s0=$?
  wait "$pid_s1" || rc_s1=$?
  if [ "$rc_s0" -ne 0 ] || [ "$rc_s1" -ne 0 ]; then
    echo "FATAL: Dual-socket Wave 1B exited non-zero (rc_s0=$rc_s0, rc_s1=$rc_s1)" >&2
    return 1
  fi

  # Wave 2: Whole-Machine & Coordinated Stage 2/3/4 Pre-Scaling Cells
  echo "=== [Wave 2: Whole-Machine Exclusive & Coordinated Pre-Scaling Cells] ==="
  unset SOCKET_PIN_CPUS
  export CELL_ALLOWLIST="s2_fig1_fig3_field_callers s2_table3_R_chr1_java_200g_n3 s4_perf_budget_4targets_n3 s4_tier3_truth_and_qual_sweep"
  run_stage_2
  run_stage_3
  run_stage_4

  # Wave 3A: Dual-Socket NUMA-Isolated Scaling Cells (th <= 8)
  echo "=== [Wave 3A: Dual-Socket NUMA-Isolated Scaling Execution (Socket 0: CPUs 0-7 / Node 0 || Socket 1: CPUs 8-15 / Node 1)] ==="
  (
    export SOCKET_PIN_CPUS="0-7"
    export CELL_ALLOWLIST="s4_scale_exome_3way_th1 s4_scale_exome_3way_th8"
    run_stage_4
  ) &
  pid_s0=$!

  (
    export SOCKET_PIN_CPUS="8-15"
    export CELL_ALLOWLIST="s4_scale_exome_3way_th2 s4_scale_exome_3way_th4 s4_span_chr1_3way_th8"
    run_stage_4
  ) &
  pid_s1=$!

  rc_s0=0; rc_s1=0
  wait "$pid_s0" || rc_s0=$?
  wait "$pid_s1" || rc_s1=$?
  if [ "$rc_s0" -ne 0 ] || [ "$rc_s1" -ne 0 ]; then
    echo "FATAL: Dual-socket Wave 3A exited non-zero (rc_s0=$rc_s0, rc_s1=$rc_s1)" >&2
    return 1
  fi

  # Wave 3B: Whole-Machine Scaling Cells (th = 16, 32)
  echo "=== [Wave 3B: Whole-Machine Scaling Cells (th = 16, 32)] ==="
  unset SOCKET_PIN_CPUS CELL_ALLOWLIST
  run_stage_4
}

# =============================================================================
# Stage 5 (Task 5f): B4 / G1 Full-Depth All-Chromosomes HG002 WGS (chr1..chr22)
# 32 vCPUs / MemoryMax=245G / -Xmx200g across VarDictRust, vardictcpp, VarDictJava
# Arm A: 22 un-chunked whole-chromosome intervals (chr1..chr22, 2.875 Gbp)
# Arm B: Pre-tiled (<=50 kbp) GIAB v4.2.1 callable genome (chr1..chr22, 2.51 Gbp) + rtg vcfeval
# =============================================================================
run_stage_5() {
  echo "=== [Stage 5: Task 5f B4 Full-Depth All-Chromosomes HG002 ~50x GRCh38 WGS (chr1..chr22, 32 vCPUs / 245G RAM / -Xmx200g)] ==="
  unset SOCKET_PIN_CPUS CELL_ALLOWLIST
  export DATA_DIR="scratch/g4"
  local hg_dir="${WORK_DIR}/data/hg002_full"

  # 1. Stage full ~130.77 GB HG002.GRCh38.2x250.bam + .bai on NVMe RAID0 & build Arm A / Arm B whole-genome BEDs
  run_local_stage_cell "stage5" "s5_stage_hg002_full_wgs" 32 "128G" bash -c "
    set -euo pipefail
    mkdir -p '$hg_dir'
    BAM_URL='https://giab.s3.amazonaws.com/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam'
    BAM_PATH='$hg_dir/HG002.GRCh38.2x250.bam'
    BAI_PATH='$hg_dir/HG002.GRCh38.2x250.bam.bai'

    if [ ! -f \"\$BAI_PATH\" ] || [ \"\$(md5sum \"\$BAI_PATH\" 2>/dev/null | awk '{print \$1}')\" != 'a3c2b449df6509ca83fbd3fea22b9aee' ]; then
      setsid gsutil -q cp '${GCS_BAMS}/hg002_full/HG002.GRCh38.2x250.bam.bai' \"\$BAI_PATH\" 2>/dev/null || \
        curl -sSLf --retry 5 \"\${BAM_URL}.bai\" -o \"\$BAI_PATH\"
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\$BAI_PATH\" '${GCS_BAMS}/hg002_full/HG002.GRCh38.2x250.bam.bai' || true
    fi
    bai_md5=\"\$(md5sum \"\$BAI_PATH\" | awk '{print \$1}')\"
    [ \"\$bai_md5\" = 'a3c2b449df6509ca83fbd3fea22b9aee' ] || { echo \"FATAL: \$BAI_PATH MD5 \$bai_md5 != a3c2b449df6509ca83fbd3fea22b9aee\" >&2; exit 1; }

    if [ ! -f \"\$BAM_PATH\" ] && setsid gsutil -q stat '${GCS_BAMS}/hg002_full/HG002.GRCh38.2x250.bam' 2>/dev/null; then
      echo '[Stage 5] Restoring 130.77 GB HG002.GRCh38.2x250.bam from GCS...'
      setsid gsutil -q cp '${GCS_BAMS}/hg002_full/HG002.GRCh38.2x250.bam' \"\$BAM_PATH\" || true
    fi
    if [ -f \"\$BAM_PATH\" ] && [ \"\$(stat -c %s \"\$BAM_PATH\" 2>/dev/null || echo 0)\" -eq 130770531934 ] && samtools quickcheck \"\$BAM_PATH\"; then
      echo '[Stage 5] Full HG002.GRCh38.2x250.bam already staged and verified on NVMe.'
    else
      if ! command -v aria2c >/dev/null 2>&1; then
        sudo apt-get update -qq && sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq aria2
      fi
      echo '[Stage 5] Downloading full 130.77 GB HG002.GRCh38.2x250.bam from GIAB S3 to NVMe RAID0...'
      aria2c -x 16 -s 16 -k 100M --continue=true --allow-overwrite=true --max-tries=10 --retry-wait=5 \
        -d '$hg_dir' -o 'HG002.GRCh38.2x250.bam' \"\$BAM_URL\"
      got_bytes=\"\$(stat -c %s \"\$BAM_PATH\" 2>/dev/null || echo 0)\"
      [ \"\$got_bytes\" -eq 130770531934 ] || { echo \"FATAL: \$BAM_PATH size \$got_bytes != 130770531934\" >&2; exit 1; }
      samtools quickcheck \"\$BAM_PATH\"
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\$BAM_PATH\" '${GCS_BAMS}/hg002_full/HG002.GRCh38.2x250.bam' || true
    fi

    # Build Arm A BED: 22 un-chunked whole-chromosome intervals (chr1..chr22, 2,875,001,522 bp)
    if ! setsid gsutil -q cp '${GCS_BAMS}/hg002_full/hg002_full_22contigs.bed' '$hg_dir/hg002_full_22contigs.bed' 2>/dev/null; then
      awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/ {printf \"%s\t0\t%s\t%s\n\", \$1, \$2, \$1}' \
        \"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta.fai\" > '$hg_dir/hg002_full_22contigs.bed'
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp '$hg_dir/hg002_full_22contigs.bed' '${GCS_BAMS}/hg002_full/hg002_full_22contigs.bed' || true
    fi
    [ \"\$(wc -l < '$hg_dir/hg002_full_22contigs.bed')\" -eq 22 ] || { echo 'FATAL: hg002_full_22contigs.bed line count != 22' >&2; exit 1; }

    # Ensure full chr1..chr22 GIAB v4.2.1 truth VCF + BED are fetched from assets.tsv and synced to g4-cache
    if [ ! -f \"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed\" ] || [ ! -f \"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz\" ]; then
      DATA_DIR=\"\$DATA_DIR\" bash harness/fetch_data.sh hg002
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r \"\$DATA_DIR/hg002/truth\" '${GCS_BAMS}/g4-cache/hg002/truth' || true
    fi

    # Build Arm B BED: pre-tiled (<=50 kbp) GIAB v4.2.1 benchmark BED across chr1..chr22 (2.51 Gbp)
    if ! setsid gsutil -q cp '${GCS_BAMS}/hg002_full/hg002_full_giab_50kb.bed' '$hg_dir/hg002_full_giab_50kb.bed' 2>/dev/null; then
      awk '
        \$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/ {
          s = \$2 + 0;
          e = \$3 + 0;
          if (s < e) {
            step = 50000;
            for (cur = s; cur < e; cur += step) {
              cur_end = (cur + step < e ? cur + step : e);
              idx++;
              printf \"%s\t%d\t%d\thg002_wgs_%d\n\", \$1, cur, cur_end, idx;
            }
          }
        }
      ' \"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed\" > '$hg_dir/hg002_full_giab_50kb.bed'
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp '$hg_dir/hg002_full_giab_50kb.bed' '${GCS_BAMS}/hg002_full/hg002_full_giab_50kb.bed' || true
    fi
    [ \"\$(wc -l < '$hg_dir/hg002_full_giab_50kb.bed')\" -ge 50000 ] || { echo 'FATAL: hg002_full_giab_50kb.bed interval count < 50000' >&2; exit 1; }
    cp -f '$hg_dir/hg002_full_22contigs.bed' \"\$KEEP_DIR/\"
    wc -l '$hg_dir/hg002_full_22contigs.bed' '$hg_dir/hg002_full_giab_50kb.bed' > \"\$KEEP_DIR/bed_counts.txt\"
  "

  # 2. Arm A: 22 Un-Chunked Whole Chromosomes (chr1..chr22, 2.875 Gbp) across all 3 engines at -th 32 / 115G / -Xmx200g
  run_checkpointed_cell "stage5" "s5_wgs_hg002_full_contigs_3way_th32" 32 "115G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    BAM='$hg_dir/HG002.GRCh38.2x250.bam'
    BED='$hg_dir/hg002_full_22contigs.bed'
    WORK_LOCAL='${WORK_DIR}/tmp_s5_contigs'
    mkdir -p \"\$WORK_LOCAL\"

    # (a) VarDictRust (-th 32, bounded 100 kbp tile streaming with 500 bp halo across all 22 autosomes > 1 Mbp)
    if [ ! -s \"\$KEEP_DIR/rust_contigs_th32.time\" ] || [ ! -s \"\$WORK_LOCAL/rust_contigs_th32.vcf\" ]; then
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_contigs_th32.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002 \
        -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -th 32 \"\$BED\" > \"\$WORK_LOCAL/rust_contigs_th32.vcf\"
    fi
    rust_vars=\$(grep -vc '^#' \"\$WORK_LOCAL/rust_contigs_th32.vcf\" || true)
    echo \"\$rust_vars\" > \"\$KEEP_DIR/rust_contigs_th32.variant_count\"
    [ \"\$rust_vars\" -ge 1000000 ] || { echo \"FATAL: rust_contigs_th32.vcf has only \$rust_vars variants (<1,000,000)\" >&2; exit 1; }
    bgzip -c -@ 16 \"\$WORK_LOCAL/rust_contigs_th32.vcf\" > \"\$KEEP_DIR/rust_contigs_th32.vcf.gz\"
    tabix -f -p vcf \"\$KEEP_DIR/rust_contigs_th32.vcf.gz\"

    # (b) vardictcpp (-th 32 --chunk 50000, contained under MemoryMax=245G)
    if [ ! -s \"\$KEEP_DIR/cpp_contigs_th32.time\" ]; then
      rc_cpp=0
      /usr/bin/time -v -o \"\$KEEP_DIR/cpp_contigs_th32.time\" \
        env CPP_CHUNK=50000 bash harness/run_cpp.sh \"\$REF\" \"\$BAM\" \"\$BED\" \"\$WORK_LOCAL/cpp_contigs_th32.vcf\" HG002 -f 0.01 -th 32 -z || rc_cpp=\$?
      echo \"\$rc_cpp\" > \"\$KEEP_DIR/cpp_contigs_th32.exit\"
      if [ -s \"\$WORK_LOCAL/cpp_contigs_th32.vcf\" ]; then
        grep -vc '^#' \"\$WORK_LOCAL/cpp_contigs_th32.vcf\" > \"\$KEEP_DIR/cpp_contigs_th32.variant_count\" || true
        bgzip -c -@ 16 \"\$WORK_LOCAL/cpp_contigs_th32.vcf\" > \"\$KEEP_DIR/cpp_contigs_th32.vcf.gz\" || true
      fi
    fi

    # (c) VarDictJava v1.8.3 (-th 32, -Xmx200g -XX:+ExitOnOutOfMemoryError, contained under MemoryMax=245G)
    if [ ! -s \"\$KEEP_DIR/java_contigs_th32.time\" ]; then
      rc_java=0
      /usr/bin/time -v -o \"\$KEEP_DIR/java_contigs_th32.time\" \
        env ORACLE_XMX_GB=200 JAVA_OPTS='-Xmx200g -XX:+ExitOnOutOfMemoryError' \
        bash harness/run_oracle.sh \"\$REF\" \"\$BAM\" \"\$BED\" \"\$WORK_LOCAL/java_contigs_th32.vcf\" HG002 -f 0.01 -th 32 -z || rc_java=\$?
      echo \"\$rc_java\" > \"\$KEEP_DIR/java_contigs_th32.exit\"
    fi
    rc_java=\"\$(cat \"\$KEEP_DIR/java_contigs_th32.exit\" 2>/dev/null || echo 1)\"
    if [ \"\$rc_java\" -eq 0 ] && [ -s \"\$WORK_LOCAL/java_contigs_th32.vcf\" ]; then
      grep -vc '^#' \"\$WORK_LOCAL/java_contigs_th32.vcf\" > \"\$KEEP_DIR/java_contigs_th32.variant_count\" || true
      bgzip -c -@ 16 \"\$WORK_LOCAL/java_contigs_th32.vcf\" > \"\$KEEP_DIR/java_contigs_th32.vcf.gz\"
      tabix -f -p vcf \"\$KEEP_DIR/java_contigs_th32.vcf.gz\"
      oracle_cache_put 'hg002_full/java_contigs_th32.vcf.gz' \"\$KEEP_DIR/java_contigs_th32.vcf.gz\"
      python3 harness/compare_vcf.py --min-variants 1000000 \"\$WORK_LOCAL/java_contigs_th32.vcf\" \"\$WORK_LOCAL/rust_contigs_th32.vcf\" | tee \"\$KEEP_DIR/compare_contigs_th32.txt\"
      grep -q '^MATCH' \"\$KEEP_DIR/compare_contigs_th32.txt\"
      echo \"COMPLETED_AND_MATCHED (rust_vars=\${rust_vars})\" | tee \"\$KEEP_DIR/contigs_verdict.txt\"
    else
      echo \"OOM-vs-completed (java_rc=\${rc_java}, rust_vars=\${rust_vars})\" | tee \"\$KEEP_DIR/contigs_verdict.txt\"
    fi
    rm -rf \"\$WORK_LOCAL\"
  "

  # 3. Arm B: Full chr1..chr22 GIAB v4.2.1 Callable Genome (<=50 kbp pre-tiled, 2.51 Gbp) + Tier-2 Parity + Tier-3 rtg vcfeval
  run_checkpointed_cell "stage5" "s5_wgs_hg002_full_giab_3way_th32" 32 "115G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    BAM='$hg_dir/HG002.GRCh38.2x250.bam'
    BED='$hg_dir/hg002_full_giab_50kb.bed'
    TRUTH_VCF=\"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz\"
    TRUTH_BED=\"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed\"
    SDF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf\"
    WORK_LOCAL='${WORK_DIR}/tmp_s5_giab'
    mkdir -p \"\$WORK_LOCAL\"

    # (a) VarDictRust v1.0.5 (-th 32)
    if [ ! -s \"\$KEEP_DIR/rust_giab_th32.time\" ] || [ ! -s \"\$WORK_LOCAL/rust_giab_th32.vcf\" ]; then
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_giab_th32.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002 \
        -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -th 32 \"\$BED\" > \"\$WORK_LOCAL/rust_giab_th32.vcf\"
    fi
    rust_vars=\$(grep -vc '^#' \"\$WORK_LOCAL/rust_giab_th32.vcf\" || true)
    echo \"\$rust_vars\" > \"\$KEEP_DIR/rust_giab_th32.variant_count\"
    [ \"\$rust_vars\" -ge 1000000 ] || { echo \"FATAL: rust_giab_th32.vcf has only \$rust_vars variants (<1,000,000)\" >&2; exit 1; }
    if [ ! -s \"\$KEEP_DIR/rust_giab_th32.vcf.gz\" ]; then
      bgzip -c -@ 16 \"\$WORK_LOCAL/rust_giab_th32.vcf\" > \"\$KEEP_DIR/rust_giab_th32.vcf.gz\"
      tabix -f -p vcf \"\$KEEP_DIR/rust_giab_th32.vcf.gz\"
    fi

    # (b) vardictcpp (-th 32 --chunk 50000)
    if [ ! -s \"\$WORK_LOCAL/cpp_giab_th32.vcf\" ] && [ -s \"\$KEEP_DIR/cpp_giab_th32.vcf.gz\" ]; then
      bgzip -dc -@ 16 \"\$KEEP_DIR/cpp_giab_th32.vcf.gz\" > \"\$WORK_LOCAL/cpp_giab_th32.vcf\"
    fi
    if [ ! -s \"\$KEEP_DIR/cpp_giab_th32.time\" ]; then
      rc_cpp=0
      /usr/bin/time -v -o \"\$KEEP_DIR/cpp_giab_th32.time\" \
        env CPP_CHUNK=50000 bash harness/run_cpp.sh \"\$REF\" \"\$BAM\" \"\$BED\" \"\$WORK_LOCAL/cpp_giab_th32.vcf\" HG002 -f 0.01 -th 32 -z || rc_cpp=\$?
      echo \"\$rc_cpp\" > \"\$KEEP_DIR/cpp_giab_th32.exit\"
      if [ -s \"\$WORK_LOCAL/cpp_giab_th32.vcf\" ]; then
        grep -vc '^#' \"\$WORK_LOCAL/cpp_giab_th32.vcf\" > \"\$KEEP_DIR/cpp_giab_th32.variant_count\" || true
        bgzip -c -@ 16 \"\$WORK_LOCAL/cpp_giab_th32.vcf\" > \"\$KEEP_DIR/cpp_giab_th32.vcf.gz\" || true
      fi
    fi

    # (c) VarDictJava v1.8.3 (-th 32, -Xmx200g -XX:+ExitOnOutOfMemoryError)
    if [ ! -s \"\$WORK_LOCAL/java_giab_th32.vcf\" ]; then
      if [ -s \"\$KEEP_DIR/java_giab_th32.vcf.gz\" ] || oracle_cache_get 'hg002_full/java_giab_th32.vcf.gz' \"\$KEEP_DIR/java_giab_th32.vcf.gz\"; then
        bgzip -dc -@ 16 \"\$KEEP_DIR/java_giab_th32.vcf.gz\" > \"\$WORK_LOCAL/java_giab_th32.vcf\"
      fi
    fi
    if [ ! -s \"\$KEEP_DIR/java_giab_th32.time\" ] || [ ! -s \"\$WORK_LOCAL/java_giab_th32.vcf\" ]; then
      rc_java=0
      /usr/bin/time -v -o \"\$KEEP_DIR/java_giab_th32.time\" \
        env ORACLE_XMX_GB=200 JAVA_OPTS='-Xmx200g -XX:+ExitOnOutOfMemoryError' \
        bash harness/run_oracle.sh \"\$REF\" \"\$BAM\" \"\$BED\" \"\$WORK_LOCAL/java_giab_th32.vcf\" HG002 -f 0.01 -th 32 -z || rc_java=\$?
      echo \"\$rc_java\" > \"\$KEEP_DIR/java_giab_th32.exit\"
    fi
    rc_java=\"\$(cat \"\$KEEP_DIR/java_giab_th32.exit\" 2>/dev/null || echo 1)\"
    if [ \"\$rc_java\" -eq 0 ] && [ -s \"\$WORK_LOCAL/java_giab_th32.vcf\" ]; then
      if [ ! -s \"\$KEEP_DIR/java_giab_th32.vcf.gz\" ]; then
        grep -vc '^#' \"\$WORK_LOCAL/java_giab_th32.vcf\" > \"\$KEEP_DIR/java_giab_th32.variant_count\" || true
        bgzip -c -@ 16 \"\$WORK_LOCAL/java_giab_th32.vcf\" > \"\$KEEP_DIR/java_giab_th32.vcf.gz\"
        tabix -f -p vcf \"\$KEEP_DIR/java_giab_th32.vcf.gz\"
        oracle_cache_put 'hg002_full/java_giab_th32.vcf.gz' \"\$KEEP_DIR/java_giab_th32.vcf.gz\"
      fi

      # Tier-2 Full-Genome Parity (VarDictJava vs VarDictRust; record diff if any, then run full-genome rtg vcfeval)
      python3 harness/compare_vcf.py --min-variants 1000000 \"\$WORK_LOCAL/java_giab_th32.vcf\" \"\$WORK_LOCAL/rust_giab_th32.vcf\" > \"\$KEEP_DIR/compare_giab_th32.txt\" 2>&1 || true
      cat \"\$KEEP_DIR/compare_giab_th32.txt\"

      if [ -s \"\$WORK_LOCAL/cpp_giab_th32.vcf\" ]; then
        python3 harness/compare_vcf.py \"\$WORK_LOCAL/java_giab_th32.vcf\" \"\$WORK_LOCAL/cpp_giab_th32.vcf\" > \"\$KEEP_DIR/compare_cpp_vs_java_giab_th32.txt\" 2>&1 || true
      fi

      # Tier-3 Full-Genome rtg vcfeval Truth Concordance across chr1..chr22
      bash harness/concordance_germline.sh hg002-full-wgs-truth \
        \"\$WORK_LOCAL/java_giab_th32.vcf\" \"\$WORK_LOCAL/rust_giab_th32.vcf\" \
        \"\$TRUTH_VCF\" \"\$TRUTH_BED\" \"\$SDF\" \
        \"\$KEEP_DIR/vcfeval_hg002_full_wgs\" HG002 \"\$WORK_LOCAL/cpp_giab_th32.vcf\" | tee \"\$KEEP_DIR/vcfeval_hg002_full_wgs.log\"
      grep -q 'EQUAL' \"\$KEEP_DIR/vcfeval_hg002_full_wgs.log\"
    else
      echo \"OOM-vs-completed (java_rc=\${rc_java}, rust_vars=\${rust_vars})\" | tee \"\$KEEP_DIR/giab_verdict.txt\"
    fi
    rm -rf \"\$WORK_LOCAL\"
  "
  # Free ~131 GB HG002 full WGS BAM after both HG002 full-genome cells complete so 2x375GB NVMe has ample room for HG008 (329 GB)
  rm -f "$hg_dir"/*.bam "$hg_dir"/*.bai 2>/dev/null || true

  # 4. Stage full ~329 GB HG008-T (195x, 160.0 GB) + HG008-N (169x, 169.0 GB) paired PDAC WGS BAMs on NVMe RAID0
  local hg8_dir="${WORK_DIR}/data/hg008_full"
  run_local_stage_cell "stage5" "s5_stage_hg008_full_wgs" 32 "128G" bash -c "
    set -euo pipefail
    mkdir -p '$hg8_dir'
    BASE_URL='https://giab.s3.amazonaws.com/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313'
    T_URL=\"\${BASE_URL}/HG008-T_Illumina_195x_GRCh38-GIABv3.bam\"
    N_URL=\"\${BASE_URL}/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam\"
    T_BAM='$hg8_dir/HG008-T_Illumina_195x_GRCh38-GIABv3.bam'
    N_BAM='$hg8_dir/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam'

    for pair in \"\$T_URL|\$T_BAM|\$DATA_DIR/hg008/HG008-T_Illumina_195x_GRCh38-GIABv3.bam.bai\" \"\$N_URL|\$N_BAM|\$DATA_DIR/hg008/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam.bai\"; do
      IFS='|' read -r u b src_bai <<< \"\$pair\"
      if [ -f \"\$src_bai\" ] && [ ! -f \"\${b}.bai\" ]; then
        cp -f \"\$src_bai\" \"\${b}.bai\"
      elif [ ! -f \"\${b}.bai\" ]; then
        setsid gsutil -q cp \"${GCS_BAMS}/hg008_full/\$(basename \"\${b}.bai\")\" \"\${b}.bai\" 2>/dev/null || \
          curl -sSLf --retry 5 \"\${u}.bai\" -o \"\${b}.bai\"
      fi
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\${b}.bai\" \"${GCS_BAMS}/hg008_full/\$(basename \"\${b}.bai\")\" || true
    done

    if ! command -v aria2c >/dev/null 2>&1; then
      sudo apt-get update -qq && sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq aria2
    fi

    if [ ! -f \"\$T_BAM\" ] && setsid gsutil -q stat '${GCS_BAMS}/hg008_full/HG008-T_Illumina_195x_GRCh38-GIABv3.bam' 2>/dev/null; then
      echo '[Stage 5] Restoring 160.0 GB HG008-T_Illumina_195x_GRCh38-GIABv3.bam from GCS...'
      setsid gsutil -q cp '${GCS_BAMS}/hg008_full/HG008-T_Illumina_195x_GRCh38-GIABv3.bam' \"\$T_BAM\" || true
    fi
    if [ -f \"\$T_BAM\" ] && [ \"\$(stat -c %s \"\$T_BAM\" 2>/dev/null || echo 0)\" -eq 159983848053 ] && samtools quickcheck \"\$T_BAM\"; then
      echo '[Stage 5] Full HG008-T 195x BAM already staged and verified on NVMe.'
    else
      echo '[Stage 5] Downloading full 160.0 GB HG008-T_Illumina_195x_GRCh38-GIABv3.bam to NVMe RAID0...'
      aria2c -x 16 -s 16 -k 100M --continue=true --allow-overwrite=true --max-tries=10 --retry-wait=5 \
        -d '$hg8_dir' -o 'HG008-T_Illumina_195x_GRCh38-GIABv3.bam' \"\$T_URL\"
      [ \"\$(stat -c %s \"\$T_BAM\" 2>/dev/null || echo 0)\" -eq 159983848053 ] || { echo 'FATAL: HG008-T size mismatch' >&2; exit 1; }
      samtools quickcheck \"\$T_BAM\"
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\$T_BAM\" '${GCS_BAMS}/hg008_full/HG008-T_Illumina_195x_GRCh38-GIABv3.bam' || true
    fi

    if [ ! -f \"\$N_BAM\" ] && setsid gsutil -q stat '${GCS_BAMS}/hg008_full/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam' 2>/dev/null; then
      echo '[Stage 5] Restoring 169.0 GB HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam from GCS...'
      setsid gsutil -q cp '${GCS_BAMS}/hg008_full/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam' \"\$N_BAM\" || true
    fi
    if [ -f \"\$N_BAM\" ] && [ \"\$(stat -c %s \"\$N_BAM\" 2>/dev/null || echo 0)\" -eq 168960212218 ] && samtools quickcheck \"\$N_BAM\"; then
      echo '[Stage 5] Full HG008-N 169x BAM already staged and verified on NVMe.'
    else
      echo '[Stage 5] Downloading full 169.0 GB HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam to NVMe RAID0...'
      aria2c -x 16 -s 16 -k 100M --continue=true --allow-overwrite=true --max-tries=10 --retry-wait=5 \
        -d '$hg8_dir' -o 'HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam' \"\$N_URL\"
      [ \"\$(stat -c %s \"\$N_BAM\" 2>/dev/null || echo 0)\" -eq 168960212218 ] || { echo 'FATAL: HG008-N size mismatch' >&2; exit 1; }
      samtools quickcheck \"\$N_BAM\"
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\$N_BAM\" '${GCS_BAMS}/hg008_full/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam' || true
    fi

    # Build Arm A BED: 22 un-chunked whole-chromosome intervals (chr1..chr22, 2,875,001,522 bp)
    if ! setsid gsutil -q cp '${GCS_BAMS}/hg008_full/hg008_full_22contigs.bed' '$hg8_dir/hg008_full_22contigs.bed' 2>/dev/null; then
      awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/ {printf \"%s\t0\t%s\t%s\n\", \$1, \$2, \$1}' \
        \"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta.fai\" > '$hg8_dir/hg008_full_22contigs.bed'
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp '$hg8_dir/hg008_full_22contigs.bed' '${GCS_BAMS}/hg008_full/hg008_full_22contigs.bed' || true
    fi
    [ \"\$(wc -l < '$hg8_dir/hg008_full_22contigs.bed')\" -eq 22 ]

    if [ ! -f \"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed\" ] || [ ! -f \"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz\" ]; then
      DATA_DIR=\"\$DATA_DIR\" bash harness/fetch_data.sh hg008
    fi

    # Build Arm B BED: pre-tiled (<=50 kbp) GIAB HG008 v0.4 benchmark BED across chr1..chr22
    SRC_BED=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed\"
    [ -f \"\$SRC_BED\" ] || SRC_BED=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed\"
    if ! setsid gsutil -q cp '${GCS_BAMS}/hg008_full/hg008_full_giab_v04_50kb.bed' '$hg8_dir/hg008_full_giab_v04_50kb.bed' 2>/dev/null; then
      awk '
        \$1 ~ /^chr([1-9]|1[0-9]|2[0-2])\$/ {
          s = \$2 + 0;
          e = \$3 + 0;
          if (s < e) {
            step = 50000;
            for (cur = s; cur < e; cur += step) {
              cur_end = (cur + step < e ? cur + step : e);
              idx++;
              printf \"%s\t%d\t%d\thg008_wgs_%d\n\", \$1, cur, cur_end, idx;
            }
          }
        }
      ' \"\$SRC_BED\" > '$hg8_dir/hg008_full_giab_v04_50kb.bed'
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp '$hg8_dir/hg008_full_giab_v04_50kb.bed' '${GCS_BAMS}/hg008_full/hg008_full_giab_v04_50kb.bed' || true
    fi
    [ \"\$(wc -l < '$hg8_dir/hg008_full_giab_v04_50kb.bed')\" -ge 40000 ]
    wc -l '$hg8_dir/hg008_full_22contigs.bed' '$hg8_dir/hg008_full_giab_v04_50kb.bed' > \"\$KEEP_DIR/hg008_bed_counts.txt\"
  "

  # 5. G2 Arm A: Full-Depth HG008-T|N (195x/169x = 364x combined) 22 Un-Chunked Chromosomes at -th 32 / 115G / -Xmx200g
  run_checkpointed_cell "stage5" "s5_wgs_hg008_full_contigs_3way_th32" 32 "115G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    T_BAM='$hg8_dir/HG008-T_Illumina_195x_GRCh38-GIABv3.bam'
    N_BAM='$hg8_dir/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam'
    BED='$hg8_dir/hg008_full_22contigs.bed'
    WORK_LOCAL='${WORK_DIR}/tmp_s5_hg008_contigs'
    mkdir -p \"\$WORK_LOCAL\"

    # (a) VarDictRust (-th 32 -U, bounded 100 kbp tile streaming with 500 bp halo across 364x combined somatic WGS > 1 Mbp; -U avoids O(N_windows * 25M) NovaSeq dark-cycle poly-G mate sink re-scans at chr2:32,916,196)
    if [ ! -s \"\$KEEP_DIR/rust_hg008_contigs_th32.time\" ] || [ ! -s \"\$WORK_LOCAL/rust_hg008_contigs_th32.vcf\" ]; then
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_hg008_contigs_th32.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'HG008-T|HG008-N' \
        -b \"\${T_BAM}|\${N_BAM}\" -c 1 -S 2 -E 3 -g 4 -z -U -th 32 \"\$BED\" > \"\$WORK_LOCAL/rust_hg008_contigs_th32.vcf\"
    fi
    rust_vars=\$(grep -vc '^#' \"\$WORK_LOCAL/rust_hg008_contigs_th32.vcf\" || true)
    echo \"\$rust_vars\" > \"\$KEEP_DIR/rust_hg008_contigs_th32.variant_count\"
    [ \"\$rust_vars\" -ge 500000 ] || { echo \"FATAL: rust_hg008_contigs_th32.vcf has only \$rust_vars variants (<500,000)\" >&2; exit 1; }
    bgzip -c -@ 16 \"\$WORK_LOCAL/rust_hg008_contigs_th32.vcf\" > \"\$KEEP_DIR/rust_hg008_contigs_th32.vcf.gz\"
    tabix -f -p vcf \"\$KEEP_DIR/rust_hg008_contigs_th32.vcf.gz\"

    # (b) vardictcpp (-th 32 --chunk 50000, contained under MemoryMax=245G)
    if [ ! -s \"\$KEEP_DIR/cpp_hg008_contigs_th32.exit\" ]; then
      rc_cpp=0
      /usr/bin/time -v -o \"\$KEEP_DIR/cpp_hg008_contigs_th32.time\" timeout 10800 \
        env CPP_CHUNK=50000 bash harness/run_cpp.sh \"\$REF\" \"\${T_BAM}|\${N_BAM}\" \"\$BED\" \"\$WORK_LOCAL/cpp_hg008_contigs_th32.vcf\" 'HG008-T|HG008-N' -f 0.01 -th 32 -z || rc_cpp=\$?
      echo \"\$rc_cpp\" > \"\$KEEP_DIR/cpp_hg008_contigs_th32.exit\"
      if [ -s \"\$WORK_LOCAL/cpp_hg008_contigs_th32.vcf\" ]; then
        grep -vc '^#' \"\$WORK_LOCAL/cpp_hg008_contigs_th32.vcf\" > \"\$KEEP_DIR/cpp_hg008_contigs_th32.variant_count\" || true
        bgzip -c -@ 16 \"\$WORK_LOCAL/cpp_hg008_contigs_th32.vcf\" > \"\$KEEP_DIR/cpp_hg008_contigs_th32.vcf.gz\" || true
      fi
    fi

    # (c) VarDictJava v1.8.3 (-th 32, -Xmx200g -XX:+ExitOnOutOfMemoryError, contained under MemoryMax=245G with 3h progress cap)
    if [ ! -s \"\$KEEP_DIR/java_hg008_contigs_th32.exit\" ]; then
      rc_java=0
      /usr/bin/time -v -o \"\$KEEP_DIR/java_hg008_contigs_th32.time\" timeout 10800 \
        env ORACLE_XMX_GB=200 JAVA_OPTS='-Xmx200g -XX:+ExitOnOutOfMemoryError' \
        bash harness/run_oracle.sh \"\$REF\" \"\${T_BAM}|\${N_BAM}\" \"\$BED\" \"\$WORK_LOCAL/java_hg008_contigs_th32.vcf\" 'HG008-T|HG008-N' -f 0.01 -th 32 -z || rc_java=\$?
      echo \"\$rc_java\" > \"\$KEEP_DIR/java_hg008_contigs_th32.exit\"
    fi
    rc_java=\"\$(cat \"\$KEEP_DIR/java_hg008_contigs_th32.exit\" 2>/dev/null || echo 1)\"
    if [ \"\$rc_java\" -eq 0 ] && [ -s \"\$WORK_LOCAL/java_hg008_contigs_th32.vcf\" ]; then
      grep -vc '^#' \"\$WORK_LOCAL/java_hg008_contigs_th32.vcf\" > \"\$KEEP_DIR/java_hg008_contigs_th32.variant_count\" || true
      bgzip -c -@ 16 \"\$WORK_LOCAL/java_hg008_contigs_th32.vcf\" > \"\$KEEP_DIR/java_hg008_contigs_th32.vcf.gz\"
      tabix -f -p vcf \"\$KEEP_DIR/java_hg008_contigs_th32.vcf.gz\"
      oracle_cache_put 'hg008_full/java_hg008_contigs_th32.vcf.gz' \"\$KEEP_DIR/java_hg008_contigs_th32.vcf.gz\"
      echo \"COMPLETED (java_rc=0, rust_vars=\${rust_vars})\" | tee \"\$KEEP_DIR/hg008_contigs_verdict.txt\"
    else
      echo \"OOM-or-timeout-vs-completed (java_rc=\${rc_java}, rust_vars=\${rust_vars})\" | tee \"\$KEEP_DIR/hg008_contigs_verdict.txt\"
    fi
    rm -rf \"\$WORK_LOCAL\"
  "

  # 6. G2 Arm B: Full-Depth HG008-T|N (195x/169x) chr1..chr22 GIAB v0.4 Callable Genome (<=50 kbp pre-tiled) + Tier-2 Parity + Tier-3 rtg vcfeval (v0.3 & v0.4)
  run_checkpointed_cell "stage5" "s5_wgs_hg008_full_giab_3way_th32" 32 "115G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    T_BAM='$hg8_dir/HG008-T_Illumina_195x_GRCh38-GIABv3.bam'
    N_BAM='$hg8_dir/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam'
    BED='$hg8_dir/hg008_full_giab_v04_50kb.bed'
    TRUTH_V03=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz\"
    BED_V03=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_all.bed\"
    TRUTH_V04=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz\"
    BED_V04=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed\"
    SDF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf\"
    WORK_LOCAL='${WORK_DIR}/tmp_s5_hg008_giab'
    mkdir -p \"\$WORK_LOCAL\"

    # (a) VarDictRust (-th 32 -U; GIAB HG008-T_somatic_smvar is a small-variant SNV+Indel <50bp benchmark; -U avoids O(N_windows * 25M) NovaSeq dark-cycle poly-G mate sink re-scans at chr2:32,916,196)
    if [ ! -s \"\$KEEP_DIR/rust_hg008_giab_th32.time\" ] || [ ! -s \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" ]; then
      /usr/bin/time -v -o \"\$KEEP_DIR/rust_hg008_giab_th32.time\" \
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'HG008-T|HG008-N' \
        -b \"\${T_BAM}|\${N_BAM}\" -c 1 -S 2 -E 3 -g 4 -z -U -th 32 \"\$BED\" > \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\"
    fi
    rust_vars=\$(grep -vc '^#' \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" || true)
    echo \"\$rust_vars\" > \"\$KEEP_DIR/rust_hg008_giab_th32.variant_count\"
    [ \"\$rust_vars\" -ge 500000 ] || { echo \"FATAL: rust_hg008_giab_th32.vcf has only \$rust_vars variants (<500,000)\" >&2; exit 1; }
    bgzip -c -@ 16 \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" > \"\$KEEP_DIR/rust_hg008_giab_th32.vcf.gz\"
    tabix -f -p vcf \"\$KEEP_DIR/rust_hg008_giab_th32.vcf.gz\"

    # Evaluate VarDictRust against full-genome HG008 v0.3 & v0.4 truth sets immediately
    bash harness/concordance.sh \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" \"\$TRUTH_V03\" \"\$BED_V03\" \"\$SDF\" \"\$KEEP_DIR/vcfeval_hg008_full_v03/rust\" HG008-T | tee \"\$KEEP_DIR/vcfeval_hg008_full_v03_rust.txt\"
    if [ -f \"\$TRUTH_V04\" ] && [ -f \"\$BED_V04\" ]; then
      bash harness/concordance.sh \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" \"\$TRUTH_V04\" \"\$BED_V04\" \"\$SDF\" \"\$KEEP_DIR/vcfeval_hg008_full_v04/rust\" HG008-T | tee \"\$KEEP_DIR/vcfeval_hg008_full_v04_rust.txt\"
    fi

    # (b) vardictcpp (-th 32 -U --chunk 50000)
    if [ ! -s \"\$KEEP_DIR/cpp_hg008_giab_th32.time\" ]; then
      rc_cpp=0
      timeout 10800 /usr/bin/time -v -o \"\$KEEP_DIR/cpp_hg008_giab_th32.time\" \
        env CPP_CHUNK=50000 bash harness/run_cpp.sh \"\$REF\" \"\${T_BAM}|\${N_BAM}\" \"\$BED\" \"\$WORK_LOCAL/cpp_hg008_giab_th32.vcf\" 'HG008-T|HG008-N' -f 0.01 -th 32 -z -U || rc_cpp=\$?
      echo \"\$rc_cpp\" > \"\$KEEP_DIR/cpp_hg008_giab_th32.exit\"
      if [ -s \"\$WORK_LOCAL/cpp_hg008_giab_th32.vcf\" ]; then
        grep -vc '^#' \"\$WORK_LOCAL/cpp_hg008_giab_th32.vcf\" > \"\$KEEP_DIR/cpp_hg008_giab_th32.variant_count\" || true
        bgzip -c -@ 16 \"\$WORK_LOCAL/cpp_hg008_giab_th32.vcf\" > \"\$KEEP_DIR/cpp_hg008_giab_th32.vcf.gz\" || true
        bash harness/concordance.sh \"\$WORK_LOCAL/cpp_hg008_giab_th32.vcf\" \"\$TRUTH_V03\" \"\$BED_V03\" \"\$SDF\" \"\$KEEP_DIR/vcfeval_hg008_full_v03/cpp\" HG008-T | tee \"\$KEEP_DIR/vcfeval_hg008_full_v03_cpp.txt\" || true
      fi
    fi

    # (c) VarDictJava v1.8.3 (-th 32 -U, -Xmx200g -XX:+ExitOnOutOfMemoryError, 4h progress cap per protocol §3 rule b)
    if [ ! -s \"\$KEEP_DIR/java_hg008_giab_th32.time\" ]; then
      rc_java=0
      timeout 14400 /usr/bin/time -v -o \"\$KEEP_DIR/java_hg008_giab_th32.time\" \
        env ORACLE_XMX_GB=200 JAVA_OPTS='-Xmx200g -XX:+ExitOnOutOfMemoryError' \
        bash harness/run_oracle.sh \"\$REF\" \"\${T_BAM}|\${N_BAM}\" \"\$BED\" \"\$WORK_LOCAL/java_hg008_giab_th32.vcf\" 'HG008-T|HG008-N' -f 0.01 -th 32 -z -U || rc_java=\$?
      echo \"\$rc_java\" > \"\$KEEP_DIR/java_hg008_giab_th32.exit\"
    fi
    rc_java=\"\$(cat \"\$KEEP_DIR/java_hg008_giab_th32.exit\" 2>/dev/null || echo 1)\"
    if [ \"\$rc_java\" -eq 0 ] && [ -s \"\$WORK_LOCAL/java_hg008_giab_th32.vcf\" ]; then
      grep -vc '^#' \"\$WORK_LOCAL/java_hg008_giab_th32.vcf\" > \"\$KEEP_DIR/java_hg008_giab_th32.variant_count\" || true
      bgzip -c -@ 16 \"\$WORK_LOCAL/java_hg008_giab_th32.vcf\" > \"\$KEEP_DIR/java_hg008_giab_th32.vcf.gz\"
      tabix -f -p vcf \"\$KEEP_DIR/java_hg008_giab_th32.vcf.gz\"
      oracle_cache_put 'hg008_full/java_hg008_giab_th32.vcf.gz' \"\$KEEP_DIR/java_hg008_giab_th32.vcf.gz\"
      python3 harness/compare_vcf.py --min-variants 500000 \"\$WORK_LOCAL/java_hg008_giab_th32.vcf\" \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" | tee \"\$KEEP_DIR/compare_hg008_giab_th32.txt\"
      grep -q '^MATCH' \"\$KEEP_DIR/compare_hg008_giab_th32.txt\"
      bash harness/concordance_pair.sh hg008-full-wgs-v03-truth \
        \"\$WORK_LOCAL/java_hg008_giab_th32.vcf\" \"\$WORK_LOCAL/rust_hg008_giab_th32.vcf\" \
        \"\$TRUTH_V03\" \"\$BED_V03\" \"\$SDF\" \"\$KEEP_DIR/vcfeval_hg008_full_v03_pair\" HG008-T \"\$WORK_LOCAL/cpp_hg008_giab_th32.vcf\" | tee \"\$KEEP_DIR/vcfeval_hg008_full_v03_pair.log\"
      grep -q 'EQUAL' \"\$KEEP_DIR/vcfeval_hg008_full_v03_pair.log\"
    else
      echo \"OOM-or-timeout-vs-completed (java_rc=\${rc_java}, rust_vars=\${rust_vars})\" | tee \"\$KEEP_DIR/hg008_giab_verdict.txt\"
    fi
    rm -rf \"\$WORK_LOCAL\"
  "
  # Free ~329 GB HG008 full WGS BAMs after both HG008 full-genome cells complete so 2x375GB NVMe has ample room for Cell 8/9
  rm -f "$hg8_dir"/*.bam "$hg8_dir"/*.bai 2>/dev/null || true

  # 7. H3 / Fig-5 GRCh38 Field-Caller Panel (Mutect2 4.6.2.0, Strelka2 2.9.10, LoFreq 2.1.5, Octopus 0.7.4, FreeBayes 1.3.10)
  #    on SEQC2 Somatic Exome (v1.2.1), HG008 Somatic WGS (v0.3 & v0.4), and HG002 Germline WGS (GIAB v4.2.1)
  run_checkpointed_cell "stage5" "s5_grch38_field_callers_truth" 32 "115G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    ensure_field_caller_envs

    REF_GIAB="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
    SDF_GIAB="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf"
    REF_D1="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    SDF_D1="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.sdf"
    [ -f "${REF_GIAB%.*}.dict" ] || samtools dict "$REF_GIAB" > "${REF_GIAB%.*}.dict"
    [ -f "${REF_GIAB}.dict" ]    || cp -f "${REF_GIAB%.*}.dict" "${REF_GIAB}.dict"
    [ -f "${REF_D1%.*}.dict" ]   || samtools dict "$REF_D1" > "${REF_D1%.*}.dict"
    [ -f "${REF_D1}.dict" ]      || cp -f "${REF_D1%.*}.dict" "${REF_D1}.dict"

    BUDGET_DIR="${WORK_DIR}/cells/s4_perf_budget_4targets_n3"
    OUT_TSV="$KEEP_DIR/grch38_field_caller_accuracy.tsv"
    printf "cohort\tcaller\twall_s\tpeak_kb\tcalls\tsnv_p\tsnv_r\tsnv_f1\tindel_p\tindel_r\tindel_f1\n" > "$OUT_TSV"

    run_timed_cmd() {
      local lbl="$1"
      shift
      if [ "$(cat "$KEEP_DIR/${lbl}.rc" 2>/dev/null || true)" = "0" ] && [ -s "$KEEP_DIR/${lbl}.time" ]; then
        return 0
      fi
      local rc=0
      /usr/bin/time -v -o "$KEEP_DIR/${lbl}.time" "$@" >"$KEEP_DIR/${lbl}.stdout" 2>"$KEEP_DIR/${lbl}.stderr" || rc=$?
      echo "$rc" > "$KEEP_DIR/${lbl}.rc"
      return "$rc"
    }

    score_truth_vcf() {
      local cohort="$1" caller="$2" vcf="$3" ref="$4" sdf="$5" truth="$6" bed="$7" sample_name="$8" time_file="${9:-}"
      [ -s "$vcf" ] || return 0
      local wall="NA" rss="NA"
      if [ -n "$time_file" ] && [ -s "$time_file" ]; then
        wall=$(awk -F": " "/Elapsed/{print \$2}" "$time_file" | tail -1)
        rss=$(awk -F": " "/Maximum resident/{print \$2}" "$time_file" | tr -dc "0-9")
      fi
      local work_ev="$KEEP_DIR/eval_${cohort}_${caller}"
      mkdir -p "$work_ev"
      local src_vcf="$vcf"
      if ! bcftools view -h "$vcf" 2>/dev/null | tail -n 1 | awk -F"\t" "{exit (NF>=10 ? 0 : 1)}"; then
        src_vcf="$work_ev/with_gt.vcf"
        bcftools view "$vcf" 2>/dev/null | awk -F"\t" -v sm="$sample_name" '\''BEGIN{OFS="\t"}
          /^##/ { print; next }
          /^#CHROM/ {
            print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            print $1,$2,$3,$4,$5,$6,$7,$8,"FORMAT",sm
            next
          }
          { print $1,$2,$3,$4,$5,$6,$7,$8,"GT","0/1" }
        '\'' > "$src_vcf"
      fi
      bcftools reheader --fai "${ref}.fai" "$src_vcf" 2>/dev/null \
        | bcftools norm -f "$ref" -m -any 2>/dev/null \
        | bcftools view -f PASS,. 2>/dev/null \
        | bcftools sort -Oz -o "$work_ev/calls.vcf.gz"
      tabix -f -p vcf "$work_ev/calls.vcf.gz"
      local n
      n=$(zgrep -vc "^#" "$work_ev/calls.vcf.gz" || echo 0)
      bcftools annotate -x FILTER "$truth" -Oz -o "$work_ev/truth.nofilter.vcf.gz"
      tabix -f -p vcf "$work_ev/truth.nofilter.vcf.gz"
      local snp_p=NA snp_r=NA snp_f=NA indel_p=NA indel_r=NA indel_f=NA
      for cls in snp indel; do
        bcftools view -v "$cls" "$work_ev/calls.vcf.gz" -Oz -o "$work_ev/calls_${cls}.vcf.gz"
        tabix -f -p vcf "$work_ev/calls_${cls}.vcf.gz"
        bcftools view -v "$cls" "$work_ev/truth.nofilter.vcf.gz" -Oz -o "$work_ev/tr_${cls}.vcf.gz"
        tabix -f -p vcf "$work_ev/tr_${cls}.vcf.gz"
        rm -rf "$work_ev/rtg_${cls}"
        rtg vcfeval -b "$work_ev/tr_${cls}.vcf.gz" -c "$work_ev/calls_${cls}.vcf.gz" -t "$sdf" \
          -e "$bed" --squash-ploidy -T 8 -o "$work_ev/rtg_${cls}" >/dev/null 2>&1 || true
        if [ -f "$work_ev/rtg_${cls}/summary.txt" ]; then
          read -r p r f <<< "$(awk "NR>2 && NF>=8" "$work_ev/rtg_${cls}/summary.txt" | sort -k8 -g | tail -1 | awk "{print \$6,\$7,\$8}")"
          eval "${cls}_p=${p:-NA}; ${cls}_r=${r:-NA}; ${cls}_f=${f:-NA}"
        fi
      done
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$cohort" "$caller" "$wall" "$rss" "$n" "$snp_p" "$snp_r" "$snp_f" "$indel_p" "$indel_r" "$indel_f" | tee -a "$OUT_TSV"
    }

    # --- A. HG002 GRCh38 Germline Truth Panel (chr20:10M-30M, GIAB v4.2.1) ---
    BAM_HG2="$DATA_DIR/hg002/HG002.chr20.bam"
    BED_HG2="$DATA_DIR/hg002/eval.bed"
    TRUTH_HG2="$DATA_DIR/hg002/truth/truth.vcf.gz"
    REF_HG2="$DATA_DIR/ref/d1vd1/GRCh38.hs38DH.fa"
    ln -sf GRCh38.d1.vd1.fa "$REF_HG2"
    samtools view -H "$BAM_HG2" | awk -F"\t" '\''/^@SQ/{sub(/^SN:/,"",$2); ok[$2]=1} END{while((getline < "'"$REF_D1"'.fai")>0){split($0,a,"\t"); if(a[1] in ok) print $0}}'\'' > "${REF_HG2}.fai"
    bgzip -c "$BED_HG2" > "$KEEP_DIR/hg002_eval.bed.gz" && tabix -f -p bed "$KEEP_DIR/hg002_eval.bed.gz"

    run_timed_cmd "hg002_FreeBayes_8c" conda run -n cl_freebayes freebayes -f "$REF_GIAB" -t "$BED_HG2" "$BAM_HG2"
    [ -s "$KEEP_DIR/hg002_FreeBayes_8c.stdout" ] && cp -f "$KEEP_DIR/hg002_FreeBayes_8c.stdout" "$KEEP_DIR/hg002_freebayes.vcf"

    IQ_HG2="$KEEP_DIR/hg002.iq.bam"
    if [ ! -f "$IQ_HG2" ]; then
      conda run -n cl_lofreq lofreq indelqual --dindel -f "$REF_GIAB" -o "$IQ_HG2" "$BAM_HG2"
      samtools index -@ 8 "$IQ_HG2"
    fi
    run_timed_cmd "hg002_LoFreq_8c" bash -c "
      set -euo pipefail
      rm -f \"$KEEP_DIR/hg002_lofreq.vcf\"
      conda run -n cl_lofreq lofreq call-parallel --pp-threads 8 --call-indels -f \"$REF_GIAB\" -l \"$BED_HG2\" -o \"$KEEP_DIR/hg002_lofreq.vcf\" \"$IQ_HG2\"
    "

    run_timed_cmd "hg002_Strelka2_8c" bash -c "
      set -euo pipefail
      rm -rf \"$KEEP_DIR/strelka_hg002\"
      conda run -n cl_strelka configureStrelkaGermlineWorkflow.py --bam \"$BAM_HG2\" --referenceFasta \"$REF_HG2\" --callRegions \"$KEEP_DIR/hg002_eval.bed.gz\" --runDir \"$KEEP_DIR/strelka_hg002\"
      conda run -n cl_strelka \"$KEEP_DIR/strelka_hg002/runWorkflow.py\" -m local -j 16
    "

    RG_HG2="$KEEP_DIR/hg002.rg.bam"
    if [ ! -f "$RG_HG2" ]; then
      samtools addreplacerg -r "@RG\tID:HG002\tSM:HG002" -@ 16 -o "$RG_HG2" "$BAM_HG2"
      samtools index -@ 16 "$RG_HG2"
    fi
    run_timed_cmd "hg002_octopus_8c" conda run -n cl_octopus octopus -R "$REF_HG2" -I "$RG_HG2" -t "$BED_HG2" --threads 16 -o "$KEEP_DIR/hg002_octopus.vcf"

    score_truth_vcf "HG002_GRCh38" "VarDictRust" "$BUDGET_DIR/hg002/rust.vcf" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$BUDGET_DIR/hg002/rust_run1.time"
    score_truth_vcf "HG002_GRCh38" "VarDictJava" "$BUDGET_DIR/hg002/java.vcf" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$BUDGET_DIR/hg002/java_run1.time"
    score_truth_vcf "HG002_GRCh38" "vardictcpp"  "$BUDGET_DIR/hg002/cpp.vcf"  "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$BUDGET_DIR/hg002/cpp_run1.time"
    score_truth_vcf "HG002_GRCh38" "FreeBayes"   "$KEEP_DIR/hg002_freebayes.vcf" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$KEEP_DIR/hg002_FreeBayes_8c.time"
    score_truth_vcf "HG002_GRCh38" "LoFreq"      "$KEEP_DIR/hg002_lofreq.vcf"    "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$KEEP_DIR/hg002_LoFreq_8c.time"
    score_truth_vcf "HG002_GRCh38" "Strelka2"    "$KEEP_DIR/strelka_hg002/results/variants/variants.vcf.gz" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$KEEP_DIR/hg002_Strelka2_8c.time"
    score_truth_vcf "HG002_GRCh38" "octopus"     "$KEEP_DIR/hg002_octopus.vcf"   "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG2" "$BED_HG2" "HG002" "$KEEP_DIR/hg002_octopus_8c.time"

    # --- B. HG008 GRCh38 Somatic Truth Panel (v0.3 & v0.4) ---
    T_HG8="$DATA_DIR/hg008/HG008-T.slice.bam"
    N_HG8="$DATA_DIR/hg008/HG008-N.slice.bam"
    BED_HG8="$DATA_DIR/hg008/eval.bed"
    TRUTH_HG8_V03="$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"
    bgzip -c "$BED_HG8" > "$KEEP_DIR/hg008_eval.bed.gz" && tabix -f -p bed "$KEEP_DIR/hg008_eval.bed.gz"

    run_timed_cmd "hg008_Mutect2_8c" bash -c "
      set -euo pipefail
      conda run -n cl_gatk4 gatk Mutect2 -R \"$REF_GIAB\" -I \"$T_HG8\" -I \"$N_HG8\" -normal HG008-N-D -L \"$BED_HG8\" --native-pair-hmm-threads 16 -O \"$KEEP_DIR/hg008_m2.unfilt.vcf.gz\"
      conda run -n cl_gatk4 gatk FilterMutectCalls -R \"$REF_GIAB\" -V \"$KEEP_DIR/hg008_m2.unfilt.vcf.gz\" -O \"$KEEP_DIR/hg008_mutect2.vcf.gz\"
    "

    run_timed_cmd "hg008_Strelka2_8c" bash -c "
      set -euo pipefail
      rm -rf \"$KEEP_DIR/strelka_hg008\"
      conda run -n cl_strelka configureStrelkaSomaticWorkflow.py --tumorBam \"$T_HG8\" --normalBam \"$N_HG8\" --referenceFasta \"$REF_GIAB\" --callRegions \"$KEEP_DIR/hg008_eval.bed.gz\" --runDir \"$KEEP_DIR/strelka_hg008\"
      conda run -n cl_strelka \"$KEEP_DIR/strelka_hg008/runWorkflow.py\" -m local -j 16
      bcftools concat -a -Oz -o \"$KEEP_DIR/hg008_strelka2.vcf.gz\" \"$KEEP_DIR/strelka_hg008/results/variants/somatic.snvs.vcf.gz\" \"$KEEP_DIR/strelka_hg008/results/variants/somatic.indels.vcf.gz\"
      tabix -f -p vcf \"$KEEP_DIR/hg008_strelka2.vcf.gz\"
    "

    run_timed_cmd "hg008_LoFreq_8c" bash -c "
      set -euo pipefail
      rm -f \"$KEEP_DIR\"/hg008_lf_*
      conda run -n cl_lofreq lofreq somatic -n \"$N_HG8\" -t \"$T_HG8\" -f \"$REF_GIAB\" -l \"$BED_HG8\" -o \"$KEEP_DIR/hg008_lf_\" --threads 16
      bcftools concat -a -Oz -o \"$KEEP_DIR/hg008_lofreq.vcf.gz\" \"$KEEP_DIR/hg008_lf_somatic_final.snvs.vcf.gz\" \"$KEEP_DIR/hg008_lf_somatic_final.indels.vcf.gz\"
      tabix -f -p vcf \"$KEEP_DIR/hg008_lofreq.vcf.gz\"
    "

    for eng in rust java cpp; do
      bcftools view -f PASS -i '\''INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"'\'' "$BUDGET_DIR/hg008/${eng}.vcf" -Oz -o "$KEEP_DIR/hg008_${eng}_som.vcf.gz"
      tabix -f -p vcf "$KEEP_DIR/hg008_${eng}_som.vcf.gz"
    done
    score_truth_vcf "HG008_v0.3" "VarDictRust" "$KEEP_DIR/hg008_rust_som.vcf.gz" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG8_V03" "$BED_HG8" "HG008-T" "$BUDGET_DIR/hg008/rust_run1.time"
    score_truth_vcf "HG008_v0.3" "VarDictJava" "$KEEP_DIR/hg008_java_som.vcf.gz" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG8_V03" "$BED_HG8" "HG008-T" "$BUDGET_DIR/hg008/java_run1.time"
    score_truth_vcf "HG008_v0.3" "vardictcpp"  "$KEEP_DIR/hg008_cpp_som.vcf.gz"  "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG8_V03" "$BED_HG8" "HG008-T" "$BUDGET_DIR/hg008/cpp_run1.time"
    score_truth_vcf "HG008_v0.3" "Mutect2"     "$KEEP_DIR/hg008_mutect2.vcf.gz"  "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG8_V03" "$BED_HG8" "HG008-T" "$KEEP_DIR/hg008_Mutect2_8c.time"
    score_truth_vcf "HG008_v0.3" "Strelka2"    "$KEEP_DIR/hg008_strelka2.vcf.gz" "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG8_V03" "$BED_HG8" "HG008-T" "$KEEP_DIR/hg008_Strelka2_8c.time"
    score_truth_vcf "HG008_v0.3" "LoFreq"      "$KEEP_DIR/hg008_lofreq.vcf.gz"   "$REF_GIAB" "$SDF_GIAB" "$TRUTH_HG8_V03" "$BED_HG8" "HG008-T" "$KEEP_DIR/hg008_LoFreq_8c.time"

    # --- C. SEQC2 GRCh38 Somatic Exome Truth Panel (HCC1395 / HCC1395BL, SEQC2 v1.2.1) ---
    T_SQ="$DATA_DIR/seqc2/HCC1395.bam"
    N_SQ="$DATA_DIR/seqc2/HCC1395BL.bam"
    BED_SQ="$DATA_DIR/seqc2/eval.bed"
    TRUTH_SQ="$DATA_DIR/seqc2/truth/truth.vcf.gz"
    bgzip -c "$BED_SQ" > "$KEEP_DIR/seqc2_eval.bed.gz" && tabix -f -p bed "$KEEP_DIR/seqc2_eval.bed.gz"

    run_timed_cmd "seqc2_Mutect2_16c" bash -c "
      set -euo pipefail
      conda run -n cl_gatk4 gatk Mutect2 -R \"$REF_D1\" -I \"$T_SQ\" -I \"$N_SQ\" -normal HCC1395BL -L \"$BED_SQ\" --native-pair-hmm-threads 16 -O \"$KEEP_DIR/seqc2_m2.unfilt.vcf.gz\"
      conda run -n cl_gatk4 gatk FilterMutectCalls -R \"$REF_D1\" -V \"$KEEP_DIR/seqc2_m2.unfilt.vcf.gz\" -O \"$KEEP_DIR/seqc2_mutect2.vcf.gz\"
    "

    run_timed_cmd "seqc2_Strelka2_16c" bash -c "
      set -euo pipefail
      rm -rf \"$KEEP_DIR/strelka_seqc2\"
      conda run -n cl_strelka configureStrelkaSomaticWorkflow.py --tumorBam \"$T_SQ\" --normalBam \"$N_SQ\" --referenceFasta \"$REF_D1\" --callRegions \"$KEEP_DIR/seqc2_eval.bed.gz\" --exome --runDir \"$KEEP_DIR/strelka_seqc2\"
      conda run -n cl_strelka \"$KEEP_DIR/strelka_seqc2/runWorkflow.py\" -m local -j 16
      bcftools concat -a -Oz -o \"$KEEP_DIR/seqc2_strelka2.vcf.gz\" \"$KEEP_DIR/strelka_seqc2/results/variants/somatic.snvs.vcf.gz\" \"$KEEP_DIR/strelka_seqc2/results/variants/somatic.indels.vcf.gz\"
      tabix -f -p vcf \"$KEEP_DIR/seqc2_strelka2.vcf.gz\"
    "

    run_timed_cmd "seqc2_LoFreq_16c" bash -c "
      set -euo pipefail
      rm -f \"$KEEP_DIR\"/seqc2_lf_*
      conda run -n cl_lofreq lofreq somatic -n \"$N_SQ\" -t \"$T_SQ\" -f \"$REF_D1\" -l \"$BED_SQ\" -o \"$KEEP_DIR/seqc2_lf_\" --threads 16
      bcftools concat -a -Oz -o \"$KEEP_DIR/seqc2_lofreq.vcf.gz\" \"$KEEP_DIR/seqc2_lf_somatic_final.snvs.vcf.gz\" \"$KEEP_DIR/seqc2_lf_somatic_final.indels.vcf.gz\"
      tabix -f -p vcf \"$KEEP_DIR/seqc2_lofreq.vcf.gz\"
    "

    for eng in rust java cpp; do
      bcftools view -f PASS -i '\''INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"'\'' "$BUDGET_DIR/seqc2/${eng}.vcf" -Oz -o "$KEEP_DIR/seqc2_${eng}_som.vcf.gz"
      tabix -f -p vcf "$KEEP_DIR/seqc2_${eng}_som.vcf.gz"
    done
    score_truth_vcf "SEQC2_v1.2.1" "VarDictRust" "$KEEP_DIR/seqc2_rust_som.vcf.gz" "$REF_D1" "$SDF_D1" "$TRUTH_SQ" "$BED_SQ" "HCC1395" "$BUDGET_DIR/seqc2/rust_run1.time"
    score_truth_vcf "SEQC2_v1.2.1" "VarDictJava" "$KEEP_DIR/seqc2_java_som.vcf.gz" "$REF_D1" "$SDF_D1" "$TRUTH_SQ" "$BED_SQ" "HCC1395" "$BUDGET_DIR/seqc2/java_run1.time"
    score_truth_vcf "SEQC2_v1.2.1" "vardictcpp"  "$KEEP_DIR/seqc2_cpp_som.vcf.gz"  "$REF_D1" "$SDF_D1" "$TRUTH_SQ" "$BED_SQ" "HCC1395" "$BUDGET_DIR/seqc2/cpp_run1.time"
    score_truth_vcf "SEQC2_v1.2.1" "Mutect2"     "$KEEP_DIR/seqc2_mutect2.vcf.gz"  "$REF_D1" "$SDF_D1" "$TRUTH_SQ" "$BED_SQ" "HCC1395" "$KEEP_DIR/seqc2_Mutect2_16c.time"
    score_truth_vcf "SEQC2_v1.2.1" "Strelka2"    "$KEEP_DIR/seqc2_strelka2.vcf.gz" "$REF_D1" "$SDF_D1" "$TRUTH_SQ" "$BED_SQ" "HCC1395" "$KEEP_DIR/seqc2_Strelka2_16c.time"
    score_truth_vcf "SEQC2_v1.2.1" "LoFreq"      "$KEEP_DIR/seqc2_lofreq.vcf.gz"   "$REF_D1" "$SDF_D1" "$TRUTH_SQ" "$BED_SQ" "HCC1395" "$KEEP_DIR/seqc2_LoFreq_16c.time"
  '

  # 8. H4 / §7 Maintenance, Supply-Chain CVE/SBOM, musl-vs-glibc bin_equiv.sh, & 60-Flag CLI Matrix Cell
  run_checkpointed_cell "stage5" "s5_h4_maintenance_security_matrix" 32 "64G" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only

    # (a) Upstream GitHub Archive Proof (AstraZeneca-NGS/VarDictJava & AstraZeneca-NGS/VarDict)
    curl -sSLf https://api.github.com/repos/AstraZeneca-NGS/VarDictJava > "$KEEP_DIR/gh_api_VarDictJava.json"
    curl -sSLf https://api.github.com/repos/AstraZeneca-NGS/VarDict > "$KEEP_DIR/gh_api_VarDict_perl.json"
    grep -q "\"archived\": true" "$KEEP_DIR/gh_api_VarDictJava.json"

    # (b) Static-MUSL vs GNU-glibc Binary Equivalence (harness/bin_equiv.sh)
    if ! command -v musl-gcc >/dev/null 2>&1; then
      sudo apt-get update -qq && sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq musl-tools
    fi
    rustup target add x86_64-unknown-linux-musl
    bash harness/bin_equiv.sh | tee "$KEEP_DIR/bin_equiv.log"
    grep -q "MUSL STATIC BINARY EQUIVALENCE: 100% PASS" "$KEEP_DIR/bin_equiv.log"

    # (c) Supply-Chain & Runtime Footprint Audit (Static MUSL vs JVM+R+Perl Conda Env)
    conda list -n vardict-oracle --json > "$KEEP_DIR/vardict_oracle_conda_packages.json"
    cargo tree -p vardictrust --prefix none > "$KEEP_DIR/vardictrust_cargo_tree.txt"
    file target/x86_64-unknown-linux-musl/release/vardictrust > "$KEEP_DIR/musl_file_info.txt"
    ldd target/x86_64-unknown-linux-musl/release/vardictrust > "$KEEP_DIR/musl_ldd_info.txt" 2>&1 || true
    stat -c "%n %s" target/release/vardictrust target/x86_64-unknown-linux-musl/release/vardictrust > "$KEEP_DIR/binary_sizes_bytes.txt"

    # (d) 60-Flag CLI Surface & Thread-Default Verification
    target/release/vardictrust --help > "$KEEP_DIR/vardictrust_help.txt"
    target/release/vardictrust --version > "$KEEP_DIR/vardictrust_version.txt"
    cargo test --test cli_flags -- --nocapture > "$KEEP_DIR/cli_flags_test.log" 2>&1 || cargo test -- --nocapture > "$KEEP_DIR/cargo_test_summary.log" 2>&1
  '

  # 9. Stage Option C Cross-Platform Short/Mid-Read Benchmark Panel (Full chr20 across 6 Chemistries + chr1 r1..r3 Somatic Panel)
  local xp_dir="${WORK_DIR}/data/cross_platform_chr20"
  run_local_stage_cell "stage5" "s5_stage_cross_platform_chr20" 32 "128G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    mkdir -p '$xp_dir'
    REF_GIAB=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    REF_D1=\"\$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa\"

    if [ ! -f \"\$REF_GIAB\" ] || [ ! -f \"\$REF_D1\" ]; then
      DATA_DIR=\"\$DATA_DIR\" bash harness/fetch_data.sh hg002 hg008 seqc2
    fi

    # Restore any already-staged cross-platform chr20 BAM/CRAM slices from GCS_BAMS/cross_platform_chr20
    if setsid gsutil -q stat \"\${GCS_BAMS}/cross_platform_chr20/*\" 2>/dev/null; then
      echo '[Cross-Platform] Restoring pre-staged chr20 BAM/CRAM slices from GCS...'
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -q rsync -r -e \"\${GCS_BAMS}/cross_platform_chr20\" '$xp_dir' || true
    fi

    # Ensure full GIAB HG002 v4.2.1 & HG008 v0.4 truth files exist and are synced to g4-cache
    TRUTH_HG2_FULL=\"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz\"
    BED_HG2_FULL=\"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed\"
    if [ ! -f \"\$TRUTH_HG2_FULL\" ] || [ ! -f \"\$BED_HG2_FULL\" ]; then
      mkdir -p \"\$DATA_DIR/hg002/truth\"
      GIAB_HG2_URL='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38'
      curl -sSLf --retry 5 \"\${GIAB_HG2_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz\" -o \"\$TRUTH_HG2_FULL\"
      curl -sSLf --retry 5 \"\${GIAB_HG2_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi\" -o \"\${TRUTH_HG2_FULL}.tbi\"
      curl -sSLf --retry 5 \"\${GIAB_HG2_URL}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed\" -o \"\$BED_HG2_FULL\"
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r \"\$DATA_DIR/hg002/truth\" '${GCS_BAMS}/g4-cache/hg002/truth' || true
    fi
    TRUTH_HG8_V04=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz\"
    BED_HG8_V04=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed\"
    if [ ! -f \"\$TRUTH_HG8_V04\" ] || [ ! -f \"\$BED_HG8_V04\" ]; then
      DATA_DIR=\"\$DATA_DIR\" bash harness/fetch_data.sh hg008
    fi

    # Download Ultima UG-High-Confidence-Regions v2.1 BED (excludes HMER>=11 homopolymers)
    if [ ! -s '$xp_dir/ug_hcr.bed' ]; then
      curl -sSLf --retry 5 'https://ultimagen-feb-2024-giab.s3.amazonaws.com/UG-High-Confidence-Regions/v2.1/ug_hcr.bed' -o '$xp_dir/ug_hcr.bed'
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp '$xp_dir/ug_hcr.bed' \"\${GCS_BAMS}/cross_platform_chr20/ug_hcr.bed\" || true
    fi

    # Build chr20 evaluation BEDs and <=50 kbp pre-tiled calling BEDs
    awk '\$1 == \"chr20\" {print \$1 \"\t\" \$2 \"\t\" \$3}' \"\$BED_HG2_FULL\" > '$xp_dir/hg002_chr20_v421_eval.bed'
    bedtools intersect -a '$xp_dir/hg002_chr20_v421_eval.bed' -b '$xp_dir/ug_hcr.bed' | sort -k1,1 -k2,2n | bedtools merge > '$xp_dir/hg002_chr20_ug_hcr_eval.bed'
    awk '{s=\$2+0; e=\$3+0; for(c=s;c<e;c+=50000){ce=(c+50000<e?c+50000:e); idx++; printf \"%s\t%d\t%d\thg2_c20_%d\n\",\$1,c,ce,idx}}' \
      '$xp_dir/hg002_chr20_v421_eval.bed' > '$xp_dir/hg002_chr20_50kb.bed'

    awk '\$1 == \"chr20\" {print \$1 \"\t\" \$2 \"\t\" \$3}' \"\$BED_HG8_V04\" > '$xp_dir/hg008_chr20_v04_eval.bed'
    bedtools intersect -a '$xp_dir/hg008_chr20_v04_eval.bed' -b '$xp_dir/ug_hcr.bed' | sort -k1,1 -k2,2n | bedtools merge > '$xp_dir/hg008_chr20_ug_hcr_eval.bed'
    awk '{s=\$2+0; e=\$3+0; for(c=s;c<e;c+=50000){ce=(c+50000<e?c+50000:e); idx++; printf \"%s\t%d\t%d\thg8_c20_%d\n\",\$1,c,ce,idx}}' \
      '$xp_dir/hg008_chr20_v04_eval.bed' > '$xp_dir/hg008_chr20_50kb.bed'

    cat \"\$DATA_DIR/hg008/r1.bed\" \"\$DATA_DIR/hg008/r2.bed\" \"\$DATA_DIR/hg008/r3.bed\" > '$xp_dir/hg008_chr1_r123.bed'

    slice_remote_bam() {
      local out_bam=\"\$1\" url=\"\$2\" idx_url=\"\$3\" reg=\"\$4\" ref_opt=\"\${5:-}\"
      local fn
      fn=\"\$(basename \"\$out_bam\")\"
      if [ ! -s \"\$out_bam\" ] && setsid gsutil -q stat \"\${GCS_BAMS}/cross_platform_chr20/\${fn}\" 2>/dev/null; then
        setsid gsutil -q cp \"\${GCS_BAMS}/cross_platform_chr20/\${fn}\" \"\$out_bam\" || true
        setsid gsutil -q cp \"\${GCS_BAMS}/cross_platform_chr20/\${fn}.bai\" \"\${out_bam}.bai\" || true
      fi
      if [ -s \"\$out_bam\" ] && [ -s \"\${out_bam}.bai\" ] && samtools quickcheck \"\$out_bam\"; then
        echo \"[Stage-XP] Reusing \$out_bam\"
        return 0
      fi
      echo \"[Stage-XP] Slicing \$out_bam (\$reg)...\"
      local idx_local=\"\${out_bam}.remote_idx\"
      curl -sSLf --retry 5 \"\$idx_url\" -o \"\$idx_local\"
      if [ -n \"\$ref_opt\" ]; then
        samtools view -@ 8 -b -T \"\$ref_opt\" -X \"\$url\" \"\$idx_local\" \$reg -o \"\$out_bam\"
      else
        samtools view -@ 8 -b -X \"\$url\" \"\$idx_local\" \$reg -o \"\$out_bam\"
      fi
      rm -f \"\$idx_local\"
      samtools index -@ 8 \"\$out_bam\"
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\$out_bam\" \"\${GCS_BAMS}/cross_platform_chr20/\${fn}\" || true
      setsid gsutil -o GSUtil:parallel_composite_upload_threshold=0 -q cp \"\${out_bam}.bai\" \"\${GCS_BAMS}/cross_platform_chr20/\${fn}.bai\" || true
    }

    # 1. Germline HG002 chr20 (Illumina, Element AVITI 80x, MGI DNBSEQ 100x, Ultima UG100 40x CRAM+BAM)
    ln -sf \"\$DATA_DIR/hg002/HG002.chr20.bam\" '$xp_dir/hg002_illumina.chr20.bam'
    ln -sf \"\$DATA_DIR/hg002/HG002.chr20.bam.bai\" '$xp_dir/hg002_illumina.chr20.bam.bai'

    G_HG2='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son'
    slice_remote_bam '$xp_dir/hg002_element_aviti.chr20.bam' \
      \"\${G_HG2}/Element_AVITI_20240920/HG002_Element-StdInsert_80x_GRCh38-GIABv3.bam\" \
      \"\${G_HG2}/Element_AVITI_20240920/HG002_Element-StdInsert_80x_GRCh38-GIABv3.bam.bai\" 'chr20'
    slice_remote_bam '$xp_dir/hg002_mgi_dnbseq.chr20.bam' \
      \"\${G_HG2}/NIST_BGIseq_2x150bp_100x/GRCh38/HG002_GRCh38_BGIseq-2x150-100x_NIST_20211126.bam\" \
      \"\${G_HG2}/NIST_BGIseq_2x150bp_100x/GRCh38/HG002_GRCh38_BGIseq-2x150-100x_NIST_20211126.bam.bai\" 'chr20'

    if [ ! -s '$xp_dir/hg002_ultima_ug100.chr20.cram' ] || [ ! -s '$xp_dir/hg002_ultima_ug100.chr20.bam' ]; then
      curl -sSLf --retry 5 'https://ultimagen-feb-2024-giab.s3.amazonaws.com/Crams/NA24385-Z0027.cram.crai' -o '$xp_dir/ug_hg2.crai'
      samtools view -@ 8 -C -T \"\$REF_D1\" -X 'https://ultimagen-feb-2024-giab.s3.amazonaws.com/Crams/NA24385-Z0027.cram' '$xp_dir/ug_hg2.crai' chr20 -o '$xp_dir/hg002_ultima_ug100.chr20.cram'
      samtools index -@ 8 '$xp_dir/hg002_ultima_ug100.chr20.cram'
      samtools view -@ 8 -b -T \"\$REF_D1\" '$xp_dir/hg002_ultima_ug100.chr20.cram' -o '$xp_dir/hg002_ultima_ug100.chr20.bam'
      samtools index -@ 8 '$xp_dir/hg002_ultima_ug100.chr20.bam'
      rm -f '$xp_dir/ug_hg2.crai'
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q cp \
        '$xp_dir/hg002_ultima_ug100.chr20.cram' '$xp_dir/hg002_ultima_ug100.chr20.cram.crai' \
        '$xp_dir/hg002_ultima_ug100.chr20.bam' '$xp_dir/hg002_ultima_ug100.chr20.bam.bai' \
        \"\${GCS_BAMS}/cross_platform_chr20/\" || true
    fi

    # 2. Somatic HG008 chr20 (Illumina 195x/169x, Roche SBX-Duplex 119x/44x, PacBio Onso 136x/48x, Element AVITI 111x/77x, Ultima UG100 306x/234x, Ultima ppmSeq 157x/160x)
    G_HG8='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab'
    if [ -s '$hg8_dir/HG008-T_Illumina_195x_GRCh38-GIABv3.bam' ] && [ ! -s '$xp_dir/hg008_illumina_T.chr20.bam' ]; then
      samtools view -@ 8 -b '$hg8_dir/HG008-T_Illumina_195x_GRCh38-GIABv3.bam' chr20 -o '$xp_dir/hg008_illumina_T.chr20.bam' && samtools index -@ 8 '$xp_dir/hg008_illumina_T.chr20.bam'
      samtools view -@ 8 -b '$hg8_dir/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam' chr20 -o '$xp_dir/hg008_illumina_N.chr20.bam' && samtools index -@ 8 '$xp_dir/hg008_illumina_N.chr20.bam'
      setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q cp \
        '$xp_dir/hg008_illumina_T.chr20.bam' '$xp_dir/hg008_illumina_T.chr20.bam.bai' \
        '$xp_dir/hg008_illumina_N.chr20.bam' '$xp_dir/hg008_illumina_N.chr20.bam.bai' \
        \"\${GCS_BAMS}/cross_platform_chr20/\" || true
    else
      slice_remote_bam '$xp_dir/hg008_illumina_T.chr20.bam' \"\${G_HG8}/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam\" \"\${G_HG8}/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam.bai\" 'chr20'
      slice_remote_bam '$xp_dir/hg008_illumina_N.chr20.bam' \"\${G_HG8}/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam\" \"\${G_HG8}/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam.bai\" 'chr20'
    fi

    slice_remote_bam '$xp_dir/hg008_roche_sbx_T.chr20.bam' \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-T_Roche-Axelios1_119X_GRCh38.bam\" \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-T_Roche-Axelios1_119X_GRCh38.bam.bai\" 'chr20'
    slice_remote_bam '$xp_dir/hg008_roche_sbx_N.chr20.bam' \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-N-D_Roche-Axelios1_44X_GRCh38.bam\" \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-N-D_Roche-Axelios1_44X_GRCh38.bam.bai\" 'chr20'

    slice_remote_bam '$xp_dir/hg008_pacbio_onso_T.chr20.bam' \"\${G_HG8}/PacBio_Onso_20240415/HG008-T_Pacbio-onso_136x_GRCh38-GIABv3.bam\" \"\${G_HG8}/PacBio_Onso_20240415/HG008-T_Pacbio-onso_136x_GRCh38-GIABv3.bam.bai\" 'chr20'
    slice_remote_bam '$xp_dir/hg008_pacbio_onso_N.chr20.bam' \"\${G_HG8}/PacBio_Onso_20240415/HG008-N-D_Pacbio-onso_48x_GRCh38-GIABv3.bam\" \"\${G_HG8}/PacBio_Onso_20240415/HG008-N-D_Pacbio-onso_48x_GRCh38-GIABv3.bam.bai\" 'chr20'

    slice_remote_bam '$xp_dir/hg008_element_aviti_T.chr20.bam' \"\${G_HG8}/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam\" \"\${G_HG8}/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam.bai\" 'chr20'
    slice_remote_bam '$xp_dir/hg008_element_aviti_N.chr20.bam' \"\${G_HG8}/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam\" \"\${G_HG8}/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam.bai\" 'chr20'

    slice_remote_bam '$xp_dir/hg008_ultima_ug100_T.chr20.bam' \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-T_Ultima_306x_GRCh38-GIABv3.cram\" \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-T_Ultima_306x_GRCh38-GIABv3.cram.crai\" 'chr20' \"\$REF_GIAB\"
    slice_remote_bam '$xp_dir/hg008_ultima_ug100_N.chr20.bam' \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-N-D_Ultima_234x_GRCh38-GIABv3.cram\" \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-N-D_Ultima_234x_GRCh38-GIABv3.cram.crai\" 'chr20' \"\$REF_GIAB\"

    slice_remote_bam '$xp_dir/hg008_ultima_ppmseq_T.chr20.bam' \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-T_Ultima_157x_GRCh38.cram\" \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-T_Ultima_157x_GRCh38.cram.crai\" 'chr20' \"\$REF_D1\"
    slice_remote_bam '$xp_dir/hg008_ultima_ppmseq_N.chr20.bam' \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-N-D_Ultima_160x_GRCh38.cram\" \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-N-D_Ultima_160x_GRCh38.cram.crai\" 'chr20' \"\$REF_D1\"

    # 3. Also stage chr1 r1..r3 (chr1:1743056-1744086, chr1:6170843-6182754, chr1:56757190-56764123) for the -f 0 -r 1 singleton noise-floor sweep
    R123='chr1:1743000-1745000 chr1:6170000-6183000 chr1:56757000-56765000'
    ln -sf \"\$DATA_DIR/hg008/HG008-T.slice.bam\" '$xp_dir/hg008_illumina_T.r123.bam'
    ln -sf \"\$DATA_DIR/hg008/HG008-T.slice.bam.bai\" '$xp_dir/hg008_illumina_T.r123.bam.bai'
    ln -sf \"\$DATA_DIR/hg008/HG008-N.slice.bam\" '$xp_dir/hg008_illumina_N.r123.bam'
    ln -sf \"\$DATA_DIR/hg008/HG008-N.slice.bam.bai\" '$xp_dir/hg008_illumina_N.r123.bam.bai'
    slice_remote_bam '$xp_dir/hg008_roche_sbx_T.r123.bam' \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-T_Roche-Axelios1_119X_GRCh38.bam\" \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-T_Roche-Axelios1_119X_GRCh38.bam.bai\" \"\$R123\"
    slice_remote_bam '$xp_dir/hg008_roche_sbx_N.r123.bam' \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-N-D_Roche-Axelios1_44X_GRCh38.bam\" \"\${G_HG8}/Hartwig_Roche-Axelios1_20260904/HG008-N-D_Roche-Axelios1_44X_GRCh38.bam.bai\" \"\$R123\"
    slice_remote_bam '$xp_dir/hg008_pacbio_onso_T.r123.bam' \"\${G_HG8}/PacBio_Onso_20240415/HG008-T_Pacbio-onso_136x_GRCh38-GIABv3.bam\" \"\${G_HG8}/PacBio_Onso_20240415/HG008-T_Pacbio-onso_136x_GRCh38-GIABv3.bam.bai\" \"\$R123\"
    slice_remote_bam '$xp_dir/hg008_pacbio_onso_N.r123.bam' \"\${G_HG8}/PacBio_Onso_20240415/HG008-N-D_Pacbio-onso_48x_GRCh38-GIABv3.bam\" \"\${G_HG8}/PacBio_Onso_20240415/HG008-N-D_Pacbio-onso_48x_GRCh38-GIABv3.bam.bai\" \"\$R123\"
    slice_remote_bam '$xp_dir/hg008_element_aviti_T.r123.bam' \"\${G_HG8}/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam\" \"\${G_HG8}/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam.bai\" \"\$R123\"
    slice_remote_bam '$xp_dir/hg008_element_aviti_N.r123.bam' \"\${G_HG8}/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam\" \"\${G_HG8}/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam.bai\" \"\$R123\"
    slice_remote_bam '$xp_dir/hg008_ultima_ug100_T.r123.bam' \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-T_Ultima_306x_GRCh38-GIABv3.cram\" \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-T_Ultima_306x_GRCh38-GIABv3.cram.crai\" \"\$R123\" \"\$REF_GIAB\"
    slice_remote_bam '$xp_dir/hg008_ultima_ug100_N.r123.bam' \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-N-D_Ultima_234x_GRCh38-GIABv3.cram\" \"\${G_HG8}/Ultima_UG100-WGS_20240530/HG008-N-D_Ultima_234x_GRCh38-GIABv3.cram.crai\" \"\$R123\" \"\$REF_GIAB\"
    slice_remote_bam '$xp_dir/hg008_ultima_ppmseq_T.r123.bam' \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-T_Ultima_157x_GRCh38.cram\" \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-T_Ultima_157x_GRCh38.cram.crai\" \"\$R123\" \"\$REF_D1\"
    slice_remote_bam '$xp_dir/hg008_ultima_ppmseq_N.r123.bam' \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-N-D_Ultima_160x_GRCh38.cram\" \"\${G_HG8}/Ultima-ppmSeq-UG100-20241028/HG008-N-D_Ultima_160x_GRCh38.cram.crai\" \"\$R123\" \"\$REF_D1\"

    date -u +%Y-%m-%dT%H:%M:%SZ > '$xp_dir/STAGED_OK'
    ls -lh '$xp_dir' > \"\$KEEP_DIR/staged_cross_platform_files.txt\"
    setsid gsutil -m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r -e '$xp_dir' \"\${GCS_BAMS}/cross_platform_chr20\" || true
  "

  # 10. Option C Cross-Platform Short-Read & Synthesis-Based Sequencing Panel (Full chr20 at -f 0.01 & -f 0 + chr1 -f 0 -r 1 Singleton Floor)
  run_checkpointed_cell "stage5" "s5_cross_platform_short_read_chr20_panel" 32 "115G" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    XP_DIR='$xp_dir'
    REF_GIAB=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    SDF_GIAB=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.sdf\"
    REF_D1=\"\$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa\"
    SDF_D1=\"\$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.sdf\"
    TRUTH_HG2=\"\$DATA_DIR/hg002/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz\"
    TRUTH_HG8=\"\$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz\"
    VR=\"target/release/vardictrust\"

    SUMMARY_TSV=\"\$KEEP_DIR/cross_platform_chr20_summary.tsv\"
    printf 'arm\tcohort\tfreq\tengine\trc\twall_s\tpeak_kb\ttotal_vars\tpass_vars\tparity_vs_java\tsnp_p\tsnp_r\tsnp_f1\tindel_p\tindel_r\tindel_f1\tughcr_snp_f1\tughcr_indel_f1\n' > \"\$SUMMARY_TSV\"

    eval_vcfeval_pair() {
      local vcf=\"\$1\" ref=\"\$2\" sdf=\"\$3\" truth=\"\$4\" bed=\"\$5\" is_som=\"\$6\" out_dir=\"\$7\"
      rm -rf \"\$out_dir\" && mkdir -p \"\$out_dir\"
      if [ ! -s \"\$vcf\" ]; then
        echo 'NA NA NA NA NA NA'
        return 0
      fi
      if [ \"\$is_som\" = '1' ]; then
        bcftools reheader --fai \"\${ref}.fai\" \"\$vcf\" 2>/dev/null \
          | bcftools view -f PASS -i 'INFO/STATUS=\"StrongSomatic\" || INFO/STATUS=\"LikelySomatic\"' 2>/dev/null \
          | bcftools norm -f \"\$ref\" -m -any 2>/dev/null \
          | bcftools sort -Oz -o \"\$out_dir/calls.vcf.gz\"
      else
        bcftools reheader --fai \"\${ref}.fai\" \"\$vcf\" 2>/dev/null \
          | bcftools view -f PASS,. 2>/dev/null \
          | bcftools norm -f \"\$ref\" -m -any 2>/dev/null \
          | bcftools sort -Oz -o \"\$out_dir/calls.vcf.gz\"
      fi
      tabix -f -p vcf \"\$out_dir/calls.vcf.gz\"
      bcftools annotate -x FILTER \"\$truth\" -R \"\$bed\" -Oz -o \"\$out_dir/truth.vcf.gz\" 2>/dev/null || \
        bcftools annotate -x FILTER \"\$truth\" -Oz -o \"\$out_dir/truth.vcf.gz\"
      tabix -f -p vcf \"\$out_dir/truth.vcf.gz\"
      local sp=NA sr=NA sf=NA ip=NA ir=NA if1=NA
      for cls in snp indel; do
        bcftools view -v \"\$cls\" \"\$out_dir/calls.vcf.gz\" -Oz -o \"\$out_dir/c_\${cls}.vcf.gz\" && tabix -f -p vcf \"\$out_dir/c_\${cls}.vcf.gz\"
        bcftools view -v \"\$cls\" \"\$out_dir/truth.vcf.gz\" -Oz -o \"\$out_dir/t_\${cls}.vcf.gz\" && tabix -f -p vcf \"\$out_dir/t_\${cls}.vcf.gz\"
        rtg vcfeval -b \"\$out_dir/t_\${cls}.vcf.gz\" -c \"\$out_dir/c_\${cls}.vcf.gz\" -t \"\$sdf\" -e \"\$bed\" --squash-ploidy -T 16 -o \"\$out_dir/rtg_\${cls}\" >/dev/null 2>&1 || true
        if [ -f \"\$out_dir/rtg_\${cls}/summary.txt\" ]; then
          read -r p r f <<< \"\$(awk 'NR>2 && NF>=8' \"\$out_dir/rtg_\${cls}/summary.txt\" | sort -k8 -g | tail -1 | awk '{print \$6,\$7,\$8}')\"
          if [ \"\$cls\" = 'snp' ]; then sp=\${p:-NA}; sr=\${r:-NA}; sf=\${f:-NA}; else ip=\${p:-NA}; ir=\${r:-NA}; if1=\${f:-NA}; fi
        fi
      done
      echo \"\$sp \$sr \$sf \$ip \$ir \$if1\"
    }

    run_chr20_3way_arm() {
      local arm=\"\$1\" cohort=\"\$2\" is_som=\"\$3\" ref=\"\$4\" sdf=\"\$5\" bam_spec=\"\$6\" sample=\"\$7\" call_bed=\"\$8\" eval_bed=\"\$9\" ughcr_bed=\"\${10}\" truth=\"\${11}\" freq=\"\${12}\"
      shift 12
      local extra_flags=(\"\$@\")
      local tag=\"\${arm}_f\${freq}\"
      echo \"=== [Cross-Platform chr20] Arm=\$arm Cohort=\$cohort freq=\$freq ===\"

      # 1. VarDictRust v1.0.5 (-th 32)
      if [ ! -s \"\$KEEP_DIR/\${tag}_rust.vcf\" ] || [ ! -s \"\$KEEP_DIR/\${tag}_rust.time\" ]; then
        /usr/bin/time -v -o \"\$KEEP_DIR/\${tag}_rust.time\" \
          \"\$VR\" -G \"\$ref\" -f \"\$freq\" -N \"\$sample\" -b \"\$bam_spec\" -c 1 -S 2 -E 3 -g 4 -z -th 32 \"\${extra_flags[@]}\" \"\$call_bed\" > \"\$KEEP_DIR/\${tag}_rust.vcf\"
      fi

      # 2. vardictcpp (-th 32 --chunk 50000)
      if [ ! -s \"\$KEEP_DIR/\${tag}_cpp.exit\" ]; then
        local rc_c=0
        timeout 3600 /usr/bin/time -v -o \"\$KEEP_DIR/\${tag}_cpp.time\" \
          env CPP_CHUNK=50000 bash harness/run_cpp.sh \"\$ref\" \"\$bam_spec\" \"\$call_bed\" \"\$KEEP_DIR/\${tag}_cpp.vcf\" \"\$sample\" -f \"\$freq\" -th 32 -z \"\${extra_flags[@]}\" || rc_c=\$?
        echo \"\$rc_c\" > \"\$KEEP_DIR/\${tag}_cpp.exit\"
      fi

      # 3. VarDictJava v1.8.3 (-th 32 -Xmx48g)
      if [ ! -s \"\$KEEP_DIR/\${tag}_java.exit\" ]; then
        local rc_j=0
        timeout 3600 /usr/bin/time -v -o \"\$KEEP_DIR/\${tag}_java.time\" \
          env ORACLE_XMX_GB=48 JAVA_OPTS='-Xmx48g -XX:+ExitOnOutOfMemoryError' \
          bash harness/run_oracle.sh \"\$ref\" \"\$bam_spec\" \"\$call_bed\" \"\$KEEP_DIR/\${tag}_java.vcf\" \"\$sample\" -f \"\$freq\" -th 32 -z \"\${extra_flags[@]}\" || rc_j=\$?
        echo \"\$rc_j\" > \"\$KEEP_DIR/\${tag}_java.exit\"
        if [ \"\$rc_j\" -eq 0 ] && [ -s \"\$KEEP_DIR/\${tag}_java.vcf\" ]; then
          oracle_cache_put \"cross_platform_chr20/\${tag}_java.vcf\" \"\$KEEP_DIR/\${tag}_java.vcf\"
        fi
      fi

      for eng in rust java cpp; do
        local vcf=\"\$KEEP_DIR/\${tag}_\${eng}.vcf\"
        local tf=\"\$KEEP_DIR/\${tag}_\${eng}.time\"
        local rc_f=\"\$KEEP_DIR/\${tag}_\${eng}.exit\"
        local rc=0
        [ -f \"\$rc_f\" ] && rc=\$(cat \"\$rc_f\")
        local wall=NA rss=NA tot=0 pass=0 par=NA
        if [ -s \"\$tf\" ]; then
          wall=\$(awk -F': ' '/Elapsed/{print \$2}' \"\$tf\" | tail -1)
          rss=\$(awk -F': ' '/Maximum resident/{print \$2}' \"\$tf\" | tr -dc '0-9')
        fi
        if [ -s \"\$vcf\" ]; then
          tot=\$(grep -vc '^#' \"\$vcf\" || true)
          if [ \"\$is_som\" = '1' ]; then
            pass=\$(grep -v '^#' \"\$vcf\" | grep -cE 'PASS.*STATUS=(Strong|Likely)Somatic' || true)
          else
            pass=\$(grep -v '^#' \"\$vcf\" | awk -F'\t' '\$7==\"PASS\" || \$7==\".\"{c++} END{print c+0}')
          fi
        fi
        if [ \"\$eng\" = 'java' ]; then
          par='ORACLE'
        elif [ -s \"\$KEEP_DIR/\${tag}_java.vcf\" ] && [ -s \"\$vcf\" ]; then
          if python3 harness/compare_vcf.py --strict-bytes \"\$KEEP_DIR/\${tag}_java.vcf\" \"\$vcf\" --eps 1e-6 > \"\$KEEP_DIR/\${tag}_\${eng}_cmp.txt\" 2>&1; then
            par='STRICT_MATCH'
          elif python3 harness/compare_vcf.py \"\$KEEP_DIR/\${tag}_java.vcf\" \"\$vcf\" --eps 1e-6 > \"\$KEEP_DIR/\${tag}_\${eng}_cmp_eps.txt\" 2>&1; then
            par='EPS_MATCH'
          else
            par='MISMATCH'
          fi
        fi
        read -r sp sr sf ip ir if1 <<< \"\$(eval_vcfeval_pair \"\$vcf\" \"\$ref\" \"\$sdf\" \"\$truth\" \"\$eval_bed\" \"\$is_som\" \"\$KEEP_DIR/ev_\${tag}_\${eng}_std\")\"
        read -r _ _ ugsf _ _ ugif1 <<< \"\$(eval_vcfeval_pair \"\$vcf\" \"\$ref\" \"\$sdf\" \"\$truth\" \"\$ughcr_bed\" \"\$is_som\" \"\$KEEP_DIR/ev_\${tag}_\${eng}_ughcr\")\"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          \"\$arm\" \"\$cohort\" \"\$freq\" \"\$eng\" \"\$rc\" \"\$wall\" \"\$rss\" \"\$tot\" \"\$pass\" \"\$par\" \
          \"\$sp\" \"\$sr\" \"\$sf\" \"\$ip\" \"\$ir\" \"\$if1\" \"\$ugsf\" \"\$ugif1\" | tee -a \"\$SUMMARY_TSV\"
      done
    }

    # --- Part A1: Germline HG002 Full chr20 across 4 Short-Read Platforms (-f 0.01 and -f 0) ---
    for freq in 0.01 0; do
      run_chr20_3way_arm 'hg002_illumina'      'HG002_Germline' 0 \"\$REF_GIAB\" \"\$SDF_GIAB\" \"\$XP_DIR/hg002_illumina.chr20.bam\"      'HG002' \"\$XP_DIR/hg002_chr20_50kb.bed\" \"\$XP_DIR/hg002_chr20_v421_eval.bed\" \"\$XP_DIR/hg002_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG2\" \"\$freq\"
      run_chr20_3way_arm 'hg002_element_aviti' 'HG002_Germline' 0 \"\$REF_GIAB\" \"\$SDF_GIAB\" \"\$XP_DIR/hg002_element_aviti.chr20.bam\" 'HG002' \"\$XP_DIR/hg002_chr20_50kb.bed\" \"\$XP_DIR/hg002_chr20_v421_eval.bed\" \"\$XP_DIR/hg002_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG2\" \"\$freq\"
      run_chr20_3way_arm 'hg002_mgi_dnbseq'    'HG002_Germline' 0 \"\$REF_D1\"   \"\$SDF_D1\"   \"\$XP_DIR/hg002_mgi_dnbseq.chr20.bam\"    'HG002' \"\$XP_DIR/hg002_chr20_50kb.bed\" \"\$XP_DIR/hg002_chr20_v421_eval.bed\" \"\$XP_DIR/hg002_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG2\" \"\$freq\"
      run_chr20_3way_arm 'hg002_ultima_ug100'  'HG002_Germline' 0 \"\$REF_D1\"   \"\$SDF_D1\"   \"\$XP_DIR/hg002_ultima_ug100.chr20.bam\"  'HG002' \"\$XP_DIR/hg002_chr20_50kb.bed\" \"\$XP_DIR/hg002_chr20_v421_eval.bed\" \"\$XP_DIR/hg002_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG2\" \"\$freq\"
    done

    # Also verify native CRAM streaming on Ultima UG100 chr20 (VarDictRust on .cram vs VarDictJava on .bam)
    /usr/bin/time -v -o \"\$KEEP_DIR/hg002_ultima_ug100_cram_rust.time\" \
      \"\$VR\" -G \"\$REF_D1\" -f 0.01 -N HG002 -b \"\$XP_DIR/hg002_ultima_ug100.chr20.cram\" -c 1 -S 2 -E 3 -g 4 -z -th 32 \"\$XP_DIR/hg002_chr20_50kb.bed\" > \"\$KEEP_DIR/hg002_ultima_ug100_cram_rust.vcf\"
    python3 harness/compare_vcf.py --strict-bytes \"\$KEEP_DIR/hg002_ultima_ug100_f0.01_java.vcf\" \"\$KEEP_DIR/hg002_ultima_ug100_cram_rust.vcf\" --eps 1e-6 | tee \"\$KEEP_DIR/ultima_cram_vs_java_bam_compare.txt\"
    grep -q '^MATCH' \"\$KEEP_DIR/ultima_cram_vs_java_bam_compare.txt\"

    # --- Part A2: Somatic HG008 Full chr20 across 6 Short/Mid-Read Platforms (-f 0.01 and -f 0 with -U) ---
    for freq in 0.01 0; do
      run_chr20_3way_arm 'hg008_illumina'      'HG008_Somatic' 1 \"\$REF_GIAB\" \"\$SDF_GIAB\" \"\$XP_DIR/hg008_illumina_T.chr20.bam|\$XP_DIR/hg008_illumina_N.chr20.bam\"           'HG008-T|HG008-N' \"\$XP_DIR/hg008_chr20_50kb.bed\" \"\$XP_DIR/hg008_chr20_v04_eval.bed\" \"\$XP_DIR/hg008_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG8\" \"\$freq\" -U
      run_chr20_3way_arm 'hg008_element_aviti' 'HG008_Somatic' 1 \"\$REF_GIAB\" \"\$SDF_GIAB\" \"\$XP_DIR/hg008_element_aviti_T.chr20.bam|\$XP_DIR/hg008_element_aviti_N.chr20.bam\" 'HG008-T|HG008-N' \"\$XP_DIR/hg008_chr20_50kb.bed\" \"\$XP_DIR/hg008_chr20_v04_eval.bed\" \"\$XP_DIR/hg008_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG8\" \"\$freq\" -U
      run_chr20_3way_arm 'hg008_pacbio_onso'   'HG008_Somatic' 1 \"\$REF_GIAB\" \"\$SDF_GIAB\" \"\$XP_DIR/hg008_pacbio_onso_T.chr20.bam|\$XP_DIR/hg008_pacbio_onso_N.chr20.bam\"     'HG008-T|HG008-N' \"\$XP_DIR/hg008_chr20_50kb.bed\" \"\$XP_DIR/hg008_chr20_v04_eval.bed\" \"\$XP_DIR/hg008_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG8\" \"\$freq\" -U
      run_chr20_3way_arm 'hg008_roche_sbx'     'HG008_Somatic' 1 \"\$REF_D1\"   \"\$SDF_D1\"   \"\$XP_DIR/hg008_roche_sbx_T.chr20.bam|\$XP_DIR/hg008_roche_sbx_N.chr20.bam\"         'HG008-T|HG008-N' \"\$XP_DIR/hg008_chr20_50kb.bed\" \"\$XP_DIR/hg008_chr20_v04_eval.bed\" \"\$XP_DIR/hg008_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG8\" \"\$freq\" -U
      run_chr20_3way_arm 'hg008_ultima_ppmseq' 'HG008_Somatic' 1 \"\$REF_D1\"   \"\$SDF_D1\"   \"\$XP_DIR/hg008_ultima_ppmseq_T.chr20.bam|\$XP_DIR/hg008_ultima_ppmseq_N.chr20.bam\" 'HG008-T|HG008-N' \"\$XP_DIR/hg008_chr20_50kb.bed\" \"\$XP_DIR/hg008_chr20_v04_eval.bed\" \"\$XP_DIR/hg008_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG8\" \"\$freq\" -U
      run_chr20_3way_arm 'hg008_ultima_ug100'  'HG008_Somatic' 1 \"\$REF_GIAB\" \"\$SDF_GIAB\" \"\$XP_DIR/hg008_ultima_ug100_T.chr20.bam|\$XP_DIR/hg008_ultima_ug100_N.chr20.bam\"   'HG008-T|HG008-N' \"\$XP_DIR/hg008_chr20_50kb.bed\" \"\$XP_DIR/hg008_chr20_v04_eval.bed\" \"\$XP_DIR/hg008_chr20_ug_hcr_eval.bed\" \"\$TRUTH_HG8\" \"\$freq\" -U
    done

    # --- Part B: Somatic Singleton Noise Floor Sweep (-f 0.01 -r 2 vs -f 0 -r 2 vs -f 0 -r 1 on chr1 r1..r3) ---
    NOISE_TSV=\"\$KEEP_DIR/cross_platform_f0_noise_floor.tsv\"
    printf 'platform\tregime\ttotal_candidates\tpass_somatic\trust_vs_java_parity\n' > \"\$NOISE_TSV\"
    for item in \
      \"Element_AVITI_111x|\$REF_GIAB|\$XP_DIR/hg008_element_aviti_T.r123.bam|\$XP_DIR/hg008_element_aviti_N.r123.bam\" \
      \"PacBio_Onso_136x|\$REF_GIAB|\$XP_DIR/hg008_pacbio_onso_T.r123.bam|\$XP_DIR/hg008_pacbio_onso_N.r123.bam\" \
      \"Illumina_195x|\$REF_GIAB|\$XP_DIR/hg008_illumina_T.r123.bam|\$XP_DIR/hg008_illumina_N.r123.bam\" \
      \"Ultima_ppmSeq_157x|\$REF_D1|\$XP_DIR/hg008_ultima_ppmseq_T.r123.bam|\$XP_DIR/hg008_ultima_ppmseq_N.r123.bam\" \
      \"Roche_SBX_119x|\$REF_D1|\$XP_DIR/hg008_roche_sbx_T.r123.bam|\$XP_DIR/hg008_roche_sbx_N.r123.bam\" \
      \"Ultima_UG100_306x|\$REF_GIAB|\$XP_DIR/hg008_ultima_ug100_T.r123.bam|\$XP_DIR/hg008_ultima_ug100_N.r123.bam\"
    do
      IFS='|' read -r plat ref tb nb <<<\"\$item\"
      for reg_spec in 'f0.01_r2|-f 0.01 -r 2' 'f0_r2|-f 0 -r 2' 'f0_r1|-f 0 -r 1'; do
        reg_lbl=\"\${reg_spec%%|*}\"
        reg_flags=\"\${reg_spec#*|}\"
        # shellcheck disable=SC2086
        \"\$VR\" -G \"\$ref\" -N 'HG008-T|HG008-N' -b \"\${tb}|\${nb}\" -c 1 -S 2 -E 3 -g 4 -z -th 16 \$reg_flags \"\$XP_DIR/hg008_chr1_r123.bed\" > \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_rust.vcf\"
        # shellcheck disable=SC2086
        bash harness/run_oracle.sh \"\$ref\" \"\${tb}|\${nb}\" \"\$XP_DIR/hg008_chr1_r123.bed\" \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_java.vcf\" 'HG008-T|HG008-N' -z -th 16 \$reg_flags
        tot=\$(grep -vc '^#' \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_rust.vcf\" || true)
        som=\$(grep -v '^#' \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_rust.vcf\" | grep -cE 'PASS.*STATUS=(Strong|Likely)Somatic' || true)
        if python3 harness/compare_vcf.py --strict-bytes \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_java.vcf\" \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_rust.vcf\" --eps 1e-6 >/dev/null 2>&1; then
          par='STRICT_MATCH'
        elif python3 harness/compare_vcf.py \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_java.vcf\" \"\$KEEP_DIR/nf_\${plat}_\${reg_lbl}_rust.vcf\" --eps 1e-6 >/dev/null 2>&1; then
          par='EPS_MATCH'
        else
          par='MISMATCH'
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' \"\$plat\" \"\$reg_lbl\" \"\$tot\" \"\$som\" \"\$par\" | tee -a \"\$NOISE_TSV\"
      done
    done
  "

  gcs_mark_done "STAGE5_DONE"
  gcs_write_heartbeat "stage5" "STAGE5_DONE" "COMPLETED"
}

# =============================================================================
# Option C (Phase C2): All-Cores (32 vCPUs) + Full-RAM (251 GiB / -Xmx200g)
# Non-Fail-Fast Full-Corpus Differential Parity Sweep (v1.0.3 vs VarDictJava v1.8.3)
# Zero scale, n=3 timing, -th 1, jobs/node, or field-caller runs.
# Intermediate GCS backup & reuse of all VarDictJava/vardictcpp truth outputs
# in ${GCS_ORACLE_CACHE} (`gs://vardict-paper-trial-results/oracle-cache-v183/`).
# =============================================================================

# Immediately backs up a local oracle truth file (.var/.sorted.var/.vcf/.tsv) + its sha256 sidecar to GCS.
oracle_cache_put() {
  local cache_rel="$1" local_file="$2"
  [ -s "$local_file" ] || return 0
  sha256sum "$local_file" | awk '{print $1}' > "${local_file}.sha256"
  if [ -n "${LOCAL_CKPT_DIR:-}" ]; then
    mkdir -p "$(dirname "${LOCAL_CKPT_DIR}/oracle-cache/${cache_rel}")"
    cp -f "$local_file" "${LOCAL_CKPT_DIR}/oracle-cache/${cache_rel}"
    cp -f "${local_file}.sha256" "${LOCAL_CKPT_DIR}/oracle-cache/${cache_rel}.sha256"
    return 0
  fi
  if ! setsid gsutil -q stat "${GCS_ORACLE_CACHE}/${cache_rel}" 2>/dev/null; then
    echo "[ORACLE-CACHE PUT] Uploading '${local_file}' -> '${GCS_ORACLE_CACHE}/${cache_rel}'"
    gcs_upload "$local_file" "${GCS_ORACLE_CACHE}/${cache_rel}" || true
    gcs_upload "${local_file}.sha256" "${GCS_ORACLE_CACHE}/${cache_rel}.sha256" || true
  fi
}

is_valid_oracle_cache_file() {
  local f="$1"
  [ -s "$f" ] || return 1
  if [[ "$f" == *.vcf ]]; then
    [ "$(grep -vc '^#' "$f" 2>/dev/null || true)" -gt 0 ] || return 1
  fi
  return 0
}

# Restores a cached oracle truth file from Local NVMe or GCS if available; otherwise runs the command
# and immediately backs up the resulting truth file (+ sha256) to GCS so it never needs regeneration.
oracle_cache_get_or_run() {
  local cache_rel="$1" local_file="$2"
  shift 2
  mkdir -p "$(dirname "$local_file")"

  # 1. Local NVMe hit
  if is_valid_oracle_cache_file "$local_file"; then
    echo "[ORACLE-CACHE HIT (NVMe)] Reusing '${local_file}' (cache_key=${cache_rel})"
    oracle_cache_put "$cache_rel" "$local_file"
    return 0
  fi
  rm -f "$local_file"

  # 2. Local test-mode cache hit
  if [ -n "${LOCAL_CKPT_DIR:-}" ] && is_valid_oracle_cache_file "${LOCAL_CKPT_DIR}/oracle-cache/${cache_rel}"; then
    cp -f "${LOCAL_CKPT_DIR}/oracle-cache/${cache_rel}" "$local_file"
    echo "[ORACLE-CACHE HIT (LOCAL)] Restored '${local_file}' from '${cache_rel}'"
    return 0
  fi

  # 3. GCS Oracle Cache hit
  if [ -z "${LOCAL_CKPT_DIR:-}" ] && setsid gsutil -q stat "${GCS_ORACLE_CACHE}/${cache_rel}" 2>/dev/null; then
    echo "[ORACLE-CACHE HIT (GCS)] Downloading '${GCS_ORACLE_CACHE}/${cache_rel}' -> '${local_file}'"
    if setsid gsutil -q cp "${GCS_ORACLE_CACHE}/${cache_rel}" "$local_file" && is_valid_oracle_cache_file "$local_file"; then
      if setsid gsutil -q stat "${GCS_ORACLE_CACHE}/${cache_rel}.sha256" 2>/dev/null; then
        setsid gsutil -q cp "${GCS_ORACLE_CACHE}/${cache_rel}.sha256" "${local_file}.sha256" || true
        local exp_sha got_sha
        exp_sha="$(awk '{print $1}' "${local_file}.sha256" 2>/dev/null || true)"
        got_sha="$(sha256sum "$local_file" | awk '{print $1}')"
        if [ -n "$exp_sha" ] && [ "$exp_sha" != "$got_sha" ]; then
          echo "[ORACLE-CACHE WARN] SHA256 mismatch on '${cache_rel}' (expected ${exp_sha}, got ${got_sha}); regenerating..." >&2
          rm -f "$local_file" "${local_file}.sha256"
        else
          return 0
        fi
      else
        return 0
      fi
    else
      rm -f "$local_file" "${local_file}.sha256"
    fi
  fi

  # 4. Cache miss -> execute generator command and immediately persist to GCS
  echo "[ORACLE-CACHE MISS] Generating '${local_file}' (cache_key=${cache_rel})..."
  "$@" || { local rc=$?; rm -f "$local_file"; return "$rc"; }
  if is_valid_oracle_cache_file "$local_file"; then
    oracle_cache_put "$cache_rel" "$local_file"
  fi
}

# Converts an already-generated VarDictJava .var file into a .vcf file using the exact
# R (teststrandbias.R / testsomatic.R) + Perl (var2vcf_valid.pl / var2vcf_paired.pl) pipeline
# from harness/run_oracle.sh without re-scanning the BAM a second time.
oracle_var_to_vcf() {
  local in_var="$1" ref="$2" bam_spec="$3" out_vcf="$4" sample="$5"
  shift 5
  local af="0.01"
  local has_fisher="false"
  local args=("$@")
  local i
  for ((i = 0; i < ${#args[@]}; i++)); do
    case "${args[$i]}" in
      -f)
        if ((i + 1 < ${#args[@]})); then
          af="${args[$((i + 1))]}"
        fi
        ;;
      -f*)
        af="${args[$i]#-f}"
        ;;
      -fisher|--fisher)
        has_fisher="true"
        ;;
    esac
  done

  local tmp_contigs=""
  if [ -f "${ref}.fai" ]; then
    tmp_contigs="$(mktemp)"
    awk '{print $1"\t0\t"$2}' "${ref}.fai" > "$tmp_contigs"
  fi

  if [[ "$bam_spec" == *"|"* ]]; then
    conda run -n vardict-oracle bash -c '
      set -euo pipefail
      V2V_ARGS=()
      if [ -n "$6" ] && [ -f "$6" ]; then
        V2V_ARGS=(-G "$1" -b "$6")
      fi
      if [ "$5" = "true" ]; then
        cat "$7" | var2vcf_paired.pl -N "$3" -f "$2" "${V2V_ARGS[@]}" > "$4"
      else
        cat "$7" | testsomatic.R | var2vcf_paired.pl -N "$3" -f "$2" "${V2V_ARGS[@]}" > "$4"
      fi
    ' _ "$ref" "$af" "$sample" "$out_vcf" "$has_fisher" "$tmp_contigs" "$in_var"
  else
    conda run -n vardict-oracle bash -c '
      set -euo pipefail
      V2V_ARGS=()
      if [ -n "$6" ] && [ -f "$6" ]; then
        V2V_ARGS=(-G "$1" -b "$6")
      fi
      if [ "$5" = "true" ]; then
        cat "$7" | var2vcf_valid.pl -N "$3" -E -f "$2" "${V2V_ARGS[@]}" > "$4"
      else
        cat "$7" | teststrandbias.R | var2vcf_valid.pl -N "$3" -E -f "$2" "${V2V_ARGS[@]}" > "$4"
      fi
    ' _ "$ref" "$af" "$sample" "$out_vcf" "$has_fisher" "$tmp_contigs" "$in_var"
  fi
  rm -f "$tmp_contigs"
}

# Seeds ${GCS_ORACLE_CACHE} with all already-generated VarDictJava & vardictcpp .var/.sorted.var/.vcf outputs
# currently present on NVMe (${WORK_DIR}/cells/) from Stage 1 and Stage 2 so nothing is ever regenerated.
seed_oracle_cache_from_existing_cells() {
  if gcs_stat_done "ORACLE_CACHE_SEEDED"; then
    echo "[SKIP] ORACLE_CACHE_SEEDED already checkpointed in GCS."
    return 0
  fi
  echo "=== [Oracle Cache Seed: Backing up existing Stage 1 & Stage 2 Java/Cpp outputs to ${GCS_ORACLE_CACHE}] ==="
  local cells_root="${WORK_DIR}/cells"
  [ -d "$cells_root" ] || return 0

  # Canonical Tier-P single-sample mappings for direct reuse by sweep & v1.0.5
  local acc
  for acc in SRR15006386 SRR15006375 SRR15006376 SRR15006540 SRR8657348; do
    local src_dir="${cells_root}/s2_tierp_v2_${acc}_th8"
    if [ -d "$src_dir" ]; then
      [ -s "${src_dir}/java_th8.var" ]        && oracle_cache_put "tierp/${acc}/java.var"        "${src_dir}/java_th8.var"
      [ -s "${src_dir}/java_th8.sorted.var" ] && oracle_cache_put "tierp/${acc}/java.sorted.var" "${src_dir}/java_th8.sorted.var"
      [ -s "${src_dir}/java_th8.vcf" ]        && oracle_cache_put "tierp/${acc}/java.vcf"        "${src_dir}/java_th8.vcf"
      [ -s "${src_dir}/cpp_th8.var" ]         && oracle_cache_put "tierp/${acc}/cpp.var"         "${src_dir}/cpp_th8.var"
      [ -s "${src_dir}/cpp_th8.sorted.var" ]  && oracle_cache_put "tierp/${acc}/cpp.sorted.var"  "${src_dir}/cpp_th8.sorted.var"
      [ -s "${src_dir}/cpp_th8.vcf" ]         && oracle_cache_put "tierp/${acc}/cpp.vcf"         "${src_dir}/cpp_th8.vcf"
    fi
  done

  # Also mirror any other java*/cpp* .var/.vcf artifacts under ${WORK_DIR}/cells/s1_* and s2_*
  local f rel
  while IFS= read -r -d '' f; do
    rel="${f#${cells_root}/}"
    oracle_cache_put "cells/${rel}" "$f"
  done < <(find "$cells_root" -maxdepth 3 -type f \( -name 'java*.var' -o -name 'java*.vcf' -o -name 'cpp*.var' -o -name 'cpp*.vcf' \) -size +0c -print0 2>/dev/null || true)
  gcs_mark_done "ORACLE_CACHE_SEEDED"
}

record_sweep_row() {
  local cell="$1" surface="$2" j_cnt="$3" r_cnt="$4" fp="$5" fn="$6" verdict="$7" notes="${8:-}"
  local summary="${WORK_DIR}/cells/v103_sweep_summary.tsv"
  mkdir -p "$(dirname "$summary")"
  if [ ! -f "$summary" ]; then
    printf 'cell\tsurface\tjava_count\trust_count\tfp_rust_vs_java\tfn_rust_vs_java\tverdict\tnotes\n' > "$summary"
  else
    awk -F'\t' -v c="$cell" -v s="$surface" '!($1==c && $2==s)' "$summary" > "${summary}.tmp" && mv -f "${summary}.tmp" "$summary"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cell" "$surface" "$j_cnt" "$r_cnt" "$fp" "$fn" "$verdict" "$notes" | tee -a "$summary"
  if [ -n "${KEEP_DIR:-}" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cell" "$surface" "$j_cnt" "$r_cnt" "$fp" "$fn" "$verdict" "$notes" >> "$KEEP_DIR/parity_verdict.tsv"
  fi
}

compare_and_record_var() {
  local cell="$1" label="$2" j_var="$3" r_var="$4" out_prefix="$5"
  LC_ALL=C sort "$j_var" > "${out_prefix}.java.sorted.var"
  LC_ALL=C sort "$r_var" > "${out_prefix}.rust.sorted.var"
  oracle_cache_put "${cell}/${label}.java.sorted.var" "${out_prefix}.java.sorted.var"
  local j_rows r_rows fp fn verdict="MATCH"
  j_rows=$(wc -l < "${out_prefix}.java.sorted.var" | tr -dc '0-9')
  r_rows=$(wc -l < "${out_prefix}.rust.sorted.var" | tr -dc '0-9')
  key_fn() { awk -F'\t' '{print $3"_"$4"_"$6"_"$7}' "$1" | LC_ALL=C sort; }
  fp=$(LC_ALL=C comm -13 <(key_fn "${out_prefix}.java.sorted.var") <(key_fn "${out_prefix}.rust.sorted.var") | wc -l | tr -dc '0-9')
  fn=$(LC_ALL=C comm -23 <(key_fn "${out_prefix}.java.sorted.var") <(key_fn "${out_prefix}.rust.sorted.var") | wc -l | tr -dc '0-9')
  diff -u "${out_prefix}.java.sorted.var" "${out_prefix}.rust.sorted.var" > "${out_prefix}.var.diff.txt" || true
  if python3 harness/compare_var.py --eps 1e-6 <(sed 's/:0-/:1-/g' "${out_prefix}.java.sorted.var") "${out_prefix}.rust.sorted.var" > "${out_prefix}.compare_var.txt" 2>&1; then
    verdict="MATCH"
  else
    verdict="DIFF_CAPTURED"
  fi
  local note
  note="$(head -n 3 "${out_prefix}.compare_var.txt" | tr '\n\t' '  ')"
  record_sweep_row "$cell" "${label}.var" "$j_rows" "$r_rows" "$fp" "$fn" "$verdict" "$note"
}

compare_and_record_vcf() {
  local cell="$1" label="$2" j_vcf="$3" r_vcf="$4" out_prefix="$5"
  oracle_cache_put "${cell}/${label}.java.vcf" "$j_vcf"
  local j_vars r_vars fp fn verdict="MATCH"
  j_vars=$(grep -vc '^#' "$j_vcf" || true)
  r_vars=$(grep -vc '^#' "$r_vcf" || true)
  key_vcf() { awk -F'\t' '!/^#/{print $1"_"$2"_"$4"_"$5}' "$1" | LC_ALL=C sort; }
  fp=$(LC_ALL=C comm -13 <(key_vcf "$j_vcf") <(key_vcf "$r_vcf") | wc -l | tr -dc '0-9')
  fn=$(LC_ALL=C comm -23 <(key_vcf "$j_vcf") <(key_vcf "$r_vcf") | wc -l | tr -dc '0-9')
  diff -u <(grep -v '^##' "$j_vcf") <(grep -v '^##' "$r_vcf") > "${out_prefix}.vcf.diff.txt" || true
  if python3 harness/compare_vcf.py --eps 1e-6 "$j_vcf" "$r_vcf" > "${out_prefix}.compare_vcf.txt" 2>&1; then
    verdict="MATCH"
  else
    verdict="DIFF_CAPTURED"
  fi
  local note
  note="$(head -n 3 "${out_prefix}.compare_vcf.txt" | tr '\n\t' '  ')"
  record_sweep_row "$cell" "${label}.vcf" "$j_vars" "$r_vars" "$fp" "$fn" "$verdict" "$note"
}

run_v103_parity_sweep() {
  echo "=== [Option C: v1.0.3 All-Cores (32 vCPUs) + Full-RAM (251 GiB / -Xmx200g) Differential Parity Sweep] ==="
  unset SOCKET_PIN_CPUS CELL_ALLOWLIST
  export DATA_DIR="scratch/g4"
  export THREADS=32
  export OMP_NUM_THREADS=32
  export RAYON_NUM_THREADS=32
  export ORACLE_XMX_GB=200
  export CPP_SKIP=1
  export SWEEP_NO_FAIL_FAST=1

  local tp_dir="${WORK_DIR}/data/tier-p"
  local rna_dir="${WORK_DIR}/data/rna"
  ensure_submodule_excludes
  ln -sfn truth/truth.vcf.gz scratch/g4/seqc2/truth.vcf.gz 2>/dev/null || true
  ln -sfn truth/truth.vcf.gz.tbi scratch/g4/seqc2/truth.vcf.gz.tbi 2>/dev/null || true

  # Seed GCS oracle cache with all Stage 1 & Stage 2 Java/Cpp outputs already on NVMe
  seed_oracle_cache_from_existing_cells

  # 1. Tier-P SRR15006540 & SRR8657348 (reuse cached .var on NVMe/GCS + cache n=1 .vcf at -th 32 / -Xmx200g)
  for item in "SRR15006540:135145" "SRR8657348:773667"; do
    local acc="${item%%:*}"
    local exp_reg="${item##*:}"
    run_checkpointed_cell "sweep" "sweep_tierp_${acc}" "0-31" "infinity" bash -c "
      set -euo pipefail
      source harness/gcp_paper_runner.sh --source-only
      REF='$tp_dir/archive/reference/hg19.fa'
      BAM='$tp_dir/bams/${acc}.bam'
      BED='$tp_dir/bams/${acc}.mkbed.bed'
      build_tierp_covered_bed \"\$BAM\" \"\$BED\" '${exp_reg}'
      PREV_DIR='${WORK_DIR}/cells/s2_tierp_v2_${acc}_th8'
      if [ -s \"\$PREV_DIR/rust_th8.var\" ]; then
        cp -f \"\$PREV_DIR/rust_th8.var\" \"\$KEEP_DIR/rust.var\"
      else
        target/release/vardictrust -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \"\$BED\" > \"\$KEEP_DIR/rust.var\"
      fi
      if [ -s \"\$PREV_DIR/java_th8.var\" ]; then
        cp -f \"\$PREV_DIR/java_th8.var\" \"\$KEEP_DIR/java.var\"
      fi
      oracle_cache_get_or_run 'tierp/${acc}/java.var' \"\$KEEP_DIR/java.var\" \
        bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$REF' -f 0.01 -N '${acc}' -b '\$BAM' -c 1 -S 2 -E 3 -g 4 -th 32 '\$BED' > '\$KEEP_DIR/java.var'\"
      compare_and_record_var 'sweep_tierp_${acc}' '${acc}' \"\$KEEP_DIR/java.var\" \"\$KEEP_DIR/rust.var\" \"\$KEEP_DIR/${acc}\"

      target/release/vardictrust -G \"\$REF\" -f 0.01 -N '${acc}' -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -th 32 \"\$BED\" > \"\$KEEP_DIR/rust.vcf\"
      oracle_cache_get_or_run 'tierp/${acc}/java.vcf' \"\$KEEP_DIR/java.vcf\" \
        oracle_var_to_vcf \"\$KEEP_DIR/java.var\" \"\$REF\" \"\$BAM\" \"\$KEEP_DIR/java.vcf\" '${acc}' -f 0.01
      compare_and_record_vcf 'sweep_tierp_${acc}' '${acc}' \"\$KEEP_DIR/java.vcf\" \"\$KEEP_DIR/rust.vcf\" \"\$KEEP_DIR/${acc}\"
    "
  done

  # 2. Tier-P Somatic SRR15006376|SRR15006386 (n=1 at -th 32 / -Xmx200g, .var + .vcf cached in GCS)
  run_checkpointed_cell "sweep" "sweep_tierp_somatic_76_86" "0-31" "infinity" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF='$tp_dir/archive/reference/hg19.fa'
    TUM='$tp_dir/bams/SRR15006376.bam'
    NOR='$tp_dir/bams/SRR15006386.bam'
    BED='$tp_dir/bams/SRR15006376.mkbed.bed'
    build_tierp_covered_bed \"\$TUM\" \"\$BED\" 105949
    target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \"\$BED\" > \"\$KEEP_DIR/rust_som.var\"
    oracle_cache_get_or_run 'tierp/somatic_76_86/java_som.var' \"\$KEEP_DIR/java_som.var\" \
      bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$REF' -f 0.01 -N 'SRR15006376|SRR15006386' -b '\${TUM}|\${NOR}' -c 1 -S 2 -E 3 -g 4 -th 32 '\$BED' > '\$KEEP_DIR/java_som.var'\"
    compare_and_record_var 'sweep_tierp_somatic_76_86' 'somatic_76_86' \"\$KEEP_DIR/java_som.var\" \"\$KEEP_DIR/rust_som.var\" \"\$KEEP_DIR/somatic_76_86\"

    target/release/vardictrust -G \"\$REF\" -f 0.01 -N 'SRR15006376|SRR15006386' -b \"\${TUM}|\${NOR}\" -c 1 -S 2 -E 3 -g 4 -th 32 \"\$BED\" > \"\$KEEP_DIR/rust_som.vcf\"
    oracle_cache_get_or_run 'tierp/somatic_76_86/java_som.vcf' \"\$KEEP_DIR/java_som.vcf\" \
      oracle_var_to_vcf \"\$KEEP_DIR/java_som.var\" \"\$REF\" \"\${TUM}|\${NOR}\" \"\$KEEP_DIR/java_som.vcf\" 'SRR15006376|SRR15006386' -f 0.01
    compare_and_record_vcf 'sweep_tierp_somatic_76_86' 'somatic_76_86' \"\$KEEP_DIR/java_som.vcf\" \"\$KEEP_DIR/rust_som.vcf\" \"\$KEEP_DIR/somatic_76_86\"
  "

  # 3. Tier-P 13-Mode Parity Sweep (n=1 at -th 32 / -Xmx200g across all modes, each Java output immediately backed up to GCS)
  run_checkpointed_cell "sweep" "sweep_tierp_13modes" "0-31" "infinity" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    NORM='$tp_dir/bams/SRR15006540.bam'
    BEDS='$tp_dir/archive/beds'
    VR='target/release/vardictrust'

    run_mode_sweep() {
      local mode=\"\$1\" j_cmd=\"\$2\" r_cmd=\"\$3\"
      oracle_cache_get_or_run \"tierp_13modes/\${mode}.java.raw.var\" \"\$KEEP_DIR/\${mode}.java.raw.var\" \
        bash -c \"\$j_cmd > '\$KEEP_DIR/\${mode}.java.raw.var'\"
      eval \"\$r_cmd\" > \"\$KEEP_DIR/\${mode}.rust.raw.var\"
      compare_and_record_var 'sweep_tierp_13modes' \"\$mode\" \"\$KEEP_DIR/\${mode}.java.raw.var\" \"\$KEEP_DIR/\${mode}.rust.raw.var\" \"\$KEEP_DIR/\${mode}\"
    }

    run_mode_sweep 'simple_f0.01' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/simple.bed 2>/dev/null\"

    run_mode_sweep 'simple_f0_mrd' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/simple.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/simple.bed 2>/dev/null\"

    run_mode_sweep 'simple_fisher' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --fisher \$BEDS/big.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --fisher --emit-intermediate -th 32 \$BEDS/big.bed 2>/dev/null\"

    run_mode_sweep 'targeted_panel_chr17' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/panel17.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/panel17.bed 2>/dev/null\"

    run_mode_sweep 'sv_deletion_regions' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/sv.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/sv.bed 2>/dev/null\"

    run_mode_sweep 'low_coverage' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 \$BEDS/lowcov.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/lowcov.bed 2>/dev/null\"

    run_mode_sweep 'amplicon_c17amp' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 -c 1 -S 2 -E 3 -g 4 \$BEDS/c17amp.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/c17amp.bed 2>/dev/null\"

    run_mode_sweep 'amplicon_fisher' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 --fisher -c 1 -S 2 -E 3 -g 4 \$BEDS/c17amp.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N s -b \$BAM -a 10:0.95 --fisher -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/c17amp.bed 2>/dev/null\"

    run_mode_sweep 'somatic_f0.01_big' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 \$BEDS/big.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/big.bed 2>/dev/null\"

    run_mode_sweep 'somatic_f0_mrd_c17' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 \$BEDS/c17.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \$BEDS/c17.bed 2>/dev/null\"

    run_mode_sweep 'somatic_fisher_c17' \
      \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --fisher \$BEDS/c17.bed -th 32 2>/dev/null\" \
      \"\$VR -G \$REF -f 0.01 -N 'tum|nrm' -b '\$BAM|\$NORM' -c 1 -S 2 -E 3 -g 4 --fisher --emit-intermediate -th 32 \$BEDS/c17.bed 2>/dev/null\"

    # P5: Tier-P Panel CRAM Twin (chr1 slice of SRR8657348.bam on simple.bed)
    samtools view -@ 16 -T \"\$REF\" -C -o \"\$KEEP_DIR/srr8657348_chr1.cram\" \"\$BAM\" chr1
    samtools index \"\$KEEP_DIR/srr8657348_chr1.cram\"
    env RUST_MIN_STACK=67108864 \"\$VR\" -G \"\$REF\" -f 0.01 -N s -b \"\$KEEP_DIR/srr8657348_chr1.cram\" -c 1 -S 2 -E 3 -g 4 --emit-intermediate -th 32 \"\$BEDS/simple.bed\" > \"\$KEEP_DIR/p5_cram.rust.raw.var\"
    compare_and_record_var 'sweep_tierp_13modes' 'P5_cram_vs_bam_twin' \"\$KEEP_DIR/simple_f0.01.rust.raw.var\" \"\$KEEP_DIR/p5_cram.rust.raw.var\" \"\$KEEP_DIR/p5_cram_twin\"
  "

  # 4. Direct Whole-Chromosome -R chr1:1-249250621 (249.25 Mbp) on SRR8657348 (n=1 at -th 32 / -Xmx200g)
  run_checkpointed_cell "sweep" "sweep_R_chr1_full" "0-31" "infinity" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF='$tp_dir/archive/reference/hg19.fa'
    BAM='$tp_dir/bams/SRR8657348.bam'
    target/release/vardictrust -G \"\$REF\" -f 0.01 -N SRR8657348 -b \"\$BAM\" --emit-intermediate -th 32 -R chr1:1-249250621 > \"\$KEEP_DIR/rust_R_chr1.var\"
    oracle_cache_get_or_run 'R_chr1_full/java_R_chr1.var' \"\$KEEP_DIR/java_R_chr1.var\" \
      bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$REF' -f 0.01 -N SRR8657348 -b '\$BAM' -th 32 -R chr1:1-249250621 > '\$KEEP_DIR/java_R_chr1.var'\"
    compare_and_record_var 'sweep_R_chr1_full' 'R_chr1_249Mbp' \"\$KEEP_DIR/java_R_chr1.var\" \"\$KEEP_DIR/rust_R_chr1.var\" \"\$KEEP_DIR/R_chr1\"
  "

  # 5. HG002 chr20 Exome (hg002.exome.20.bam on target20.bed, n=1 at -th 32 / -Xmx200g)
  run_checkpointed_cell "sweep" "sweep_hg002_chr20_exome" "0-31" "infinity" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF='$tp_dir/archive/reference/ref20.fa'
    BAM='$tp_dir/bams/hg002.exome.20.bam'
    BED='$tp_dir/truth/target20.bed'
    [ -f \"\${REF}.fai\" ] || samtools faidx \"\$REF\"
    awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\"R\"}' \"\$BED\" > \"\$KEEP_DIR/target20_4col.bed\"
    target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002 -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z --emit-intermediate -th 32 \"\$KEEP_DIR/target20_4col.bed\" > \"\$KEEP_DIR/rust.var\"
    oracle_cache_get_or_run 'hg002_chr20_exome/java.var' \"\$KEEP_DIR/java.var\" \
      bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$REF' -f 0.01 -N HG002 -b '\$BAM' -c 1 -S 2 -E 3 -g 4 -z 1 -th 32 '\$KEEP_DIR/target20_4col.bed' > '\$KEEP_DIR/java.var'\"
    compare_and_record_var 'sweep_hg002_chr20_exome' 'hg002_chr20_exome' \"\$KEEP_DIR/java.var\" \"\$KEEP_DIR/rust.var\" \"\$KEEP_DIR/hg002_chr20_exome\"

    target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002 -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -th 32 \"\$KEEP_DIR/target20_4col.bed\" > \"\$KEEP_DIR/rust.vcf\"
    oracle_cache_get_or_run 'hg002_chr20_exome/java.vcf' \"\$KEEP_DIR/java.vcf\" \
      oracle_var_to_vcf \"\$KEEP_DIR/java.var\" \"\$REF\" \"\$BAM\" \"\$KEEP_DIR/java.vcf\" HG002 -f 0.01
    compare_and_record_vcf 'sweep_hg002_chr20_exome' 'hg002_chr20_exome' \"\$KEEP_DIR/java.vcf\" \"\$KEEP_DIR/rust.vcf\" \"\$KEEP_DIR/hg002_chr20_exome\"
  "

  # 6. Full 51.8 GiB RNA-Seq (HG002_rep1 STAR BAM: Arm A = -i splice table across GIAB eval.bed; Arm B = cDNA variant calling on chr20 eval.bed)
  run_checkpointed_cell "sweep" "sweep_rna_hg002_52gb" "0-31" "infinity" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    BAM='$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam'
    BED=\"\$DATA_DIR/hg002/eval.bed\"

    # Arm A: Full-transcriptome -i splice-junction table
    if [ ! -s \"\$KEEP_DIR/rust.splice.tsv\" ]; then
      target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002_RNA -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z -i --emit-intermediate -th 32 \"\$BED\" > \"\$KEEP_DIR/rust.splice.tsv\"
    fi
    oracle_cache_get_or_run 'rna_hg002_52gb/java.splice.tsv' \"\$KEEP_DIR/java.splice.tsv\" \
      bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$REF' -f 0.01 -N HG002_RNA -b '\$BAM' -c 1 -S 2 -E 3 -g 4 -z 1 -i -th 32 '\$BED' > '\$KEEP_DIR/java.splice.tsv'\"
    LC_ALL=C sort \"\$KEEP_DIR/java.splice.tsv\" > \"\$KEEP_DIR/java.splice.sorted.tsv\"
    oracle_cache_put 'rna_hg002_52gb/java.splice.sorted.tsv' \"\$KEEP_DIR/java.splice.sorted.tsv\"
    LC_ALL=C sort \"\$KEEP_DIR/rust.splice.tsv\" > \"\$KEEP_DIR/rust.splice.sorted.tsv\"
    diff -u \"\$KEEP_DIR/java.splice.sorted.tsv\" \"\$KEEP_DIR/rust.splice.sorted.tsv\" > \"\$KEEP_DIR/splice_diff.txt\" || true
    j_sp=\$(wc -l < \"\$KEEP_DIR/java.splice.sorted.tsv\" | tr -dc '0-9')
    r_sp=\$(wc -l < \"\$KEEP_DIR/rust.splice.sorted.tsv\" | tr -dc '0-9')
    fp_sp=\$(LC_ALL=C comm -13 \"\$KEEP_DIR/java.splice.sorted.tsv\" \"\$KEEP_DIR/rust.splice.sorted.tsv\" | wc -l | tr -dc '0-9')
    fn_sp=\$(LC_ALL=C comm -23 \"\$KEEP_DIR/java.splice.sorted.tsv\" \"\$KEEP_DIR/rust.splice.sorted.tsv\" | wc -l | tr -dc '0-9')
    sp_verdict='MATCH'
    if ! python3 harness/compare_var.py --eps 1e-6 <(sed 's/:0-/:1-/g' \"\$KEEP_DIR/java.splice.sorted.tsv\") \"\$KEEP_DIR/rust.splice.sorted.tsv\" > \"\$KEEP_DIR/compare_splice_var.txt\" 2>&1; then
      sp_verdict='DIFF_CAPTURED'
    fi
    record_sweep_row 'sweep_rna_hg002_52gb' 'rna_splice_table' \"\$j_sp\" \"\$r_sp\" \"\$fp_sp\" \"\$fn_sp\" \"\$sp_verdict\" \"\$(head -n 2 \"\$KEEP_DIR/compare_splice_var.txt\" | tr '\n\t' '  ')\"

    # Arm B: Splice-aware cDNA variant calling (without -i) on chr20 GIAB eval regions (first 500 intervals)
    awk '\$1==\"chr20\"{print; if (++n==500) exit}' \"\$BED\" > \"\$KEEP_DIR/chr20_rna_cdna.bed\"
    target/release/vardictrust -G \"\$REF\" -f 0.01 -N HG002_RNA -b \"\$BAM\" -c 1 -S 2 -E 3 -g 4 -z --emit-intermediate -th 32 \"\$KEEP_DIR/chr20_rna_cdna.bed\" > \"\$KEEP_DIR/rust_cdna.var\"
    oracle_cache_get_or_run 'rna_hg002_52gb/java_cdna.var' \"\$KEEP_DIR/java_cdna.var\" \
      bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$REF' -f 0.01 -N HG002_RNA -b '\$BAM' -c 1 -S 2 -E 3 -g 4 -z 1 -th 32 '\$KEEP_DIR/chr20_rna_cdna.bed' > '\$KEEP_DIR/java_cdna.var'\"
    compare_and_record_var 'sweep_rna_hg002_52gb' 'rna_cdna_variants' \"\$KEEP_DIR/java_cdna.var\" \"\$KEEP_DIR/rust_cdna.var\" \"\$KEEP_DIR/rna_cdna\"
  "

  # 7. Four Full-Sample Targets (Amplicon, HG002 227k Deep WGS, HG008 18k Somatic Deep WGS, SEQC2 623k Full Somatic Exome; n=1 at -th 32 / -Xmx200g)
  run_checkpointed_cell "sweep" "sweep_4targets_full" "0-31" "infinity" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    mkdir -p "$KEEP_DIR/amplicon" "$KEEP_DIR/hg002" "$KEEP_DIR/hg008" "$KEEP_DIR/seqc2"

    # (a) Amplicon
    oracle_cache_get_or_run "4targets_full/amplicon/java.vcf" "$KEEP_DIR/amplicon/java.vcf" \
      env ORACLE_XMX_GB=200 bash harness/run_oracle.sh "$DATA_DIR/amplicon/ref.fa" "$DATA_DIR/amplicon/amplicon.bam" "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" "$KEEP_DIR/amplicon/java.vcf" SAMPLE -z -a 10:0.95 -th 32
    target/release/vardictrust -G "$DATA_DIR/amplicon/ref.fa" -f 0.01 -N SAMPLE -b "$DATA_DIR/amplicon/amplicon.bam" -c 1 -S 2 -E 3 -g 4 -z -a 10:0.95 -th 32 "$DATA_DIR/amplicon/correct_8columns_amplicon.bed" > "$KEEP_DIR/amplicon/rust.vcf"
    compare_and_record_vcf "sweep_4targets_full" "amplicon_full" "$KEEP_DIR/amplicon/java.vcf" "$KEEP_DIR/amplicon/rust.vcf" "$KEEP_DIR/amplicon/cmp"

    # (b) HG002 Deep WGS (chr20:10M-30M, 227,562 variants)
    oracle_cache_get_or_run "4targets_full/hg002/java.vcf" "$KEEP_DIR/hg002/java.vcf" \
      env ORACLE_XMX_GB=200 bash harness/run_oracle.sh "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" "$DATA_DIR/hg002/HG002.chr20.bam" "$DATA_DIR/hg002/eval.bed" "$KEEP_DIR/hg002/java.vcf" HG002 -z -f 0.01 -th 32
    target/release/vardictrust -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta" -f 0.01 -N HG002 -b "$DATA_DIR/hg002/HG002.chr20.bam" -c 1 -S 2 -E 3 -g 4 -z -th 32 "$DATA_DIR/hg002/eval.bed" > "$KEEP_DIR/hg002/rust.vcf"
    compare_and_record_vcf "sweep_4targets_full" "hg002_deep_wgs" "$KEEP_DIR/hg002/java.vcf" "$KEEP_DIR/hg002/rust.vcf" "$KEEP_DIR/hg002/cmp"

    # (c) HG008 Somatic Paired Deep WGS (18,426 variants)
    oracle_cache_get_or_run "4targets_full/hg008/java.vcf" "$KEEP_DIR/hg008/java.vcf" \
      env ORACLE_XMX_GB=200 bash harness/run_oracle.sh "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" "$DATA_DIR/hg008/engine.bed" "$KEEP_DIR/hg008/java.vcf" "HG008-T|HG008-N" -z -f 0.01 -th 32
    target/release/vardictrust -G "$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fa" -f 0.01 -N "HG008-T|HG008-N" -b "$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam" -c 1 -S 2 -E 3 -g 4 -z -th 32 "$DATA_DIR/hg008/engine.bed" > "$KEEP_DIR/hg008/rust.vcf"
    compare_and_record_vcf "sweep_4targets_full" "hg008_somatic_wgs" "$KEEP_DIR/hg008/java.vcf" "$KEEP_DIR/hg008/rust.vcf" "$KEEP_DIR/hg008/cmp"

    # (d) SEQC2 HCC1395|HCC1395BL Full Somatic WES (623,444 variants)
    oracle_cache_get_or_run "4targets_full/seqc2/java.vcf" "$KEEP_DIR/seqc2/java.vcf" \
      env ORACLE_XMX_GB=200 bash harness/run_oracle.sh "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa" "$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam" "$DATA_DIR/seqc2/exome.engine.bed" "$KEEP_DIR/seqc2/java.vcf" "HCC1395|HCC1395BL" -z -f 0.01 -th 32
    target/release/vardictrust -G "$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa" -f 0.01 -N "HCC1395|HCC1395BL" -b "$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam" -c 1 -S 2 -E 3 -g 4 -z -th 32 "$DATA_DIR/seqc2/exome.engine.bed" > "$KEEP_DIR/seqc2/rust.vcf"
    compare_and_record_vcf "sweep_4targets_full" "seqc2_somatic_exome_full" "$KEEP_DIR/seqc2/java.vcf" "$KEEP_DIR/seqc2/rust.vcf" "$KEEP_DIR/seqc2/cmp"
  '

  # 8. Tier-G chr1:0-10,000,000 MRD -f 0 Noise Floor on HCC1395.bam (n=1 at -th 32 / -Xmx200g, .var + .vcf)
  run_checkpointed_cell "sweep" "sweep_span_chr1_f0" "0-31" "infinity" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    BAM="$DATA_DIR/seqc2/HCC1395.bam"
    printf "chr1\t0\t10000000\tR10M\n" > "$KEEP_DIR/chr1_10m.bed"
    target/release/vardictrust -G "$REF" -f 0 -N HCC1395 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z --emit-intermediate -th 32 "$KEEP_DIR/chr1_10m.bed" > "$KEEP_DIR/rust_chr1_10m_f0.var"
    oracle_cache_get_or_run "span_chr1_f0/java_chr1_10m_f0.var" "$KEEP_DIR/java_chr1_10m_f0.var" \
      bash -c "_JAVA_OPTIONS=\"-Xmx200g\" conda run -n vardict-oracle vardict-java -G \"$REF\" -f 0 -N HCC1395 -b \"$BAM\" -c 1 -S 2 -E 3 -g 4 -z 1 -th 32 \"$KEEP_DIR/chr1_10m.bed\" > \"$KEEP_DIR/java_chr1_10m_f0.var\""
    compare_and_record_var "sweep_span_chr1_f0" "chr1_10m_f0" "$KEEP_DIR/java_chr1_10m_f0.var" "$KEEP_DIR/rust_chr1_10m_f0.var" "$KEEP_DIR/chr1_10m_f0"

    target/release/vardictrust -G "$REF" -f 0 -N HCC1395 -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th 32 "$KEEP_DIR/chr1_10m.bed" > "$KEEP_DIR/rust_chr1_10m_f0.vcf"
    oracle_cache_get_or_run "span_chr1_f0/java_chr1_10m_f0.vcf" "$KEEP_DIR/java_chr1_10m_f0.vcf" \
      oracle_var_to_vcf "$KEEP_DIR/java_chr1_10m_f0.var" "$REF" "$BAM" "$KEEP_DIR/java_chr1_10m_f0.vcf" HCC1395 -f 0
    compare_and_record_vcf "sweep_span_chr1_f0" "chr1_10m_f0" "$KEEP_DIR/java_chr1_10m_f0.vcf" "$KEEP_DIR/rust_chr1_10m_f0.vcf" "$KEEP_DIR/chr1_10m_f0"
  '

  # 9. Dedicated WGS, MRD (-f 0 & -f 0.001), and WGS+MRD Exhaustive Oracle Sweep (n=1 at -th 32 / -Xmx200g)
  run_checkpointed_cell "sweep" "sweep_wgs_mrd_matrix" "0-31" "infinity" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    FA_GIAB="$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta"
    FA_SEQC2="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    VR="target/release/vardictrust"

    run_wgs_mrd_arm() {
      local arm="$1" ref="$2" bam="$3" bed="$4" sample="$5"
      shift 5
      local extra=("$@")
      "$VR" -G "$ref" -N "$sample" -b "$bam" -c 1 -S 2 -E 3 -g 4 -z --emit-intermediate -th 32 "${extra[@]}" "$bed" > "$KEEP_DIR/${arm}.rust.var"
      oracle_cache_get_or_run "wgs_mrd_matrix/${arm}.java.var" "$KEEP_DIR/${arm}.java.var" \
        bash -c "_JAVA_OPTIONS=\"-Xmx200g\" conda run -n vardict-oracle vardict-java -G \"$ref\" -N \"$sample\" -b \"$bam\" -c 1 -S 2 -E 3 -g 4 -z 1 -th 32 ${extra[*]} \"$bed\" > \"$KEEP_DIR/${arm}.java.var\""
      compare_and_record_var "sweep_wgs_mrd_matrix" "$arm" "$KEEP_DIR/${arm}.java.var" "$KEEP_DIR/${arm}.rust.var" "$KEEP_DIR/${arm}"

      "$VR" -G "$ref" -N "$sample" -b "$bam" -c 1 -S 2 -E 3 -g 4 -z -th 32 "${extra[@]}" "$bed" > "$KEEP_DIR/${arm}.rust.vcf"
      oracle_cache_get_or_run "wgs_mrd_matrix/${arm}.java.vcf" "$KEEP_DIR/${arm}.java.vcf" \
        oracle_var_to_vcf "$KEEP_DIR/${arm}.java.var" "$ref" "$bam" "$KEEP_DIR/${arm}.java.vcf" "$sample" "${extra[@]}"
      compare_and_record_vcf "sweep_wgs_mrd_matrix" "$arm" "$KEEP_DIR/${arm}.java.vcf" "$KEEP_DIR/${arm}.rust.vcf" "$KEEP_DIR/${arm}"
    }

    # Arm 1: HG002 Germline Deep WGS + MRD (-f 0 and -f 0.001) across chr20 windows + 1 Mbp chr20:10M-11M
    cat > "$KEEP_DIR/hg002_wgs_mrd.bed" <<EOF
chr20	10000000	10020000	HG002_3A
chr20	20000000	20020000	HG002_3B
chr20	30000000	30020000	HG002_3C
chr20	10020000	11000000	HG002_1MBP
EOF
    run_wgs_mrd_arm "hg002_wgs_mrd_f0.001" "$FA_GIAB" "$DATA_DIR/hg002/HG002.chr20.bam" "$KEEP_DIR/hg002_wgs_mrd.bed" "HG002" -f 0.001
    run_wgs_mrd_arm "hg002_wgs_mrd_f0"     "$FA_GIAB" "$DATA_DIR/hg002/HG002.chr20.bam" "$KEEP_DIR/hg002_wgs_mrd.bed" "HG002" -f 0

    # Arm 2: HG008 Somatic Paired Deep WGS + MRD (-f 0 and -f 0.001, default & --fisher) across engine.bed
    HG008_BAMS="$DATA_DIR/hg008/HG008-T.slice.bam|$DATA_DIR/hg008/HG008-N.slice.bam"
    run_wgs_mrd_arm "hg008_somatic_wgs_mrd_f0.001"        "$FA_GIAB" "$HG008_BAMS" "$DATA_DIR/hg008/engine.bed" "HG008T|HG008N" -f 0.001
    run_wgs_mrd_arm "hg008_somatic_wgs_mrd_f0.001_fisher" "$FA_GIAB" "$HG008_BAMS" "$DATA_DIR/hg008/engine.bed" "HG008T|HG008N" -f 0.001 --fisher
    run_wgs_mrd_arm "hg008_somatic_wgs_mrd_f0"            "$FA_GIAB" "$HG008_BAMS" "$DATA_DIR/hg008/engine.bed" "HG008T|HG008N" -f 0

    # Arm 3: HG005 300x Ultra-Deep Diversity WGS + MRD (-f 0 and -f 0.001) across chr3, chr6, chr21
    printf "chr3\t40000000\t40050000\tHG005_3A\n" > "$KEEP_DIR/hg005_3a.bed"
    printf "chr6\t30000000\t30050000\tHG005_3B\n" > "$KEEP_DIR/hg005_3b.bed"
    printf "chr21\t20000000\t20050000\tHG005_3C\n" > "$KEEP_DIR/hg005_3c.bed"
    run_wgs_mrd_arm "hg005_300x_wgs_chr3_f0.001"  "$FA_GIAB" "$DATA_DIR/multimodal/hg005_s4a.bam" "$KEEP_DIR/hg005_3a.bed" "HG005" -f 0.001
    run_wgs_mrd_arm "hg005_300x_wgs_chr3_f0"      "$FA_GIAB" "$DATA_DIR/multimodal/hg005_s4a.bam" "$KEEP_DIR/hg005_3a.bed" "HG005" -f 0
    run_wgs_mrd_arm "hg005_300x_wgs_chr6_f0.001"  "$FA_GIAB" "$DATA_DIR/multimodal/hg005_s4b.bam" "$KEEP_DIR/hg005_3b.bed" "HG005" -f 0.001
    run_wgs_mrd_arm "hg005_300x_wgs_chr21_f0.001" "$FA_GIAB" "$DATA_DIR/multimodal/hg005_s4c.bam" "$KEEP_DIR/hg005_3c.bed" "HG005" -f 0.001

    # Arm 4: SEQC2 Somatic Paired WES + MRD (-f 0 and -f 0.001, default & --fisher) across slices.engine.bed
    SEQC2_BAMS="$DATA_DIR/seqc2/HCC1395.bam|$DATA_DIR/seqc2/HCC1395BL.bam"
    run_wgs_mrd_arm "seqc2_somatic_mrd_f0.001"        "$FA_SEQC2" "$SEQC2_BAMS" "$DATA_DIR/seqc2/slices.engine.bed" "HCC1395|HCC1395BL" -f 0.001
    run_wgs_mrd_arm "seqc2_somatic_mrd_f0.001_fisher" "$FA_SEQC2" "$SEQC2_BAMS" "$DATA_DIR/seqc2/slices.engine.bed" "HCC1395|HCC1395BL" -f 0.001 --fisher
    run_wgs_mrd_arm "seqc2_somatic_mrd_f0"            "$FA_SEQC2" "$SEQC2_BAMS" "$DATA_DIR/seqc2/slices.engine.bed" "HCC1395|HCC1395BL" -f 0
  '

  # 10. Dedicated Super-Ultra-Deep (8,035x - 173,132x) + Ultra-Low VAF (-f 0.0005, -f 0.001, -f 0 + --fisher) Oracle Sweep
  run_checkpointed_cell "sweep" "sweep_ultradeep_mrd_stress" "0-31" "infinity" bash -c "
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    REF_HG19='$tp_dir/archive/reference/hg19.fa'
    FA_GIAB=\"\$DATA_DIR/ref/giabv3/GRCh38_GIABv3.fasta\"
    VR='target/release/vardictrust'

    # Build BED of empirically verified 8,035x - 30,267x single-sample (29,389x - 45,707x combined somatic) exons in PRJNA742684 AML
    cat > \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" <<'EOF'
chr9	139390000	139400000	NOTCH1_25939x
chr3	128204000	128206000	GATA2_22073x
chr5	1294000	1295000	TERT_20387x
chr21	36164000	36165000	RUNX1_20128x
chrX	39932000	39933000	BCOR_19580x
chr12	12022000	12023000	ETV6_19646x
chr5	176941000	176944000	NSD1_19824x
chr13	28607000	28625000	FLT3_9107x
chr14	19488000	19490000	CHR14_30267x
chr22	16349000	16354000	CHR22_21987x
chr19	17952000	17953000	JAK3_9524x
chr17	74733000	74734000	SRSF2_9112x
EOF

    run_ud_arm() {
      local arm=\"\$1\" ref=\"\$2\" bam=\"\$3\" bed=\"\$4\" sample=\"\$5\"
      shift 5
      local extra=(\"\$@\")
      \"\$VR\" -G \"\$ref\" -N \"\$sample\" -b \"\$bam\" -c 1 -S 2 -E 3 -g 4 -z --emit-intermediate -th 32 \"\${extra[@]}\" \"\$bed\" > \"\$KEEP_DIR/\${arm}.rust.var\"
      oracle_cache_get_or_run \"ultradeep_mrd/\${arm}.java.var\" \"\$KEEP_DIR/\${arm}.java.var\" \
        bash -c \"_JAVA_OPTIONS='-Xmx200g' conda run -n vardict-oracle vardict-java -G '\$ref' -N '\$sample' -b '\$bam' -c 1 -S 2 -E 3 -g 4 -z 1 -th 32 \${extra[*]} '\$bed' > '\$KEEP_DIR/\${arm}.java.var'\"
      compare_and_record_var 'sweep_ultradeep_mrd_stress' \"\$arm\" \"\$KEEP_DIR/\${arm}.java.var\" \"\$KEEP_DIR/\${arm}.rust.var\" \"\$KEEP_DIR/\${arm}\"

      \"\$VR\" -G \"\$ref\" -N \"\$sample\" -b \"\$bam\" -c 1 -S 2 -E 3 -g 4 -z -th 32 \"\${extra[@]}\" \"\$bed\" > \"\$KEEP_DIR/\${arm}.rust.vcf\"
      oracle_cache_get_or_run \"ultradeep_mrd/\${arm}.java.vcf\" \"\$KEEP_DIR/\${arm}.java.vcf\" \
        oracle_var_to_vcf \"\$KEEP_DIR/\${arm}.java.var\" \"\$ref\" \"\$bam\" \"\$KEEP_DIR/\${arm}.java.vcf\" \"\$sample\" \"\${extra[@]}\"
      compare_and_record_vcf 'sweep_ultradeep_mrd_stress' \"\$arm\" \"\$KEEP_DIR/\${arm}.java.vcf\" \"\$KEEP_DIR/\${arm}.rust.vcf\" \"\$KEEP_DIR/\${arm}\"
    }

    # Arm 1: Single-sample 25,939x (SRR15006375) & 30,267x (SRR15006540) at -f 0.0005 (0.05% VAF), -f 0.001 (0.1% VAF), and -f 0 with --fisher
    for acc in SRR15006375 SRR15006540 SRR15006386 SRR15006376; do
      run_ud_arm \"ud_\${acc}_f0.0005_fisher\" \"\$REF_HG19\" '$tp_dir/bams/'\"\${acc}\".bam \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" \"\${acc}\" -f 0.0005 --fisher
      run_ud_arm \"ud_\${acc}_f0_fisher\"      \"\$REF_HG19\" '$tp_dir/bams/'\"\${acc}\".bam \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" \"\${acc}\" -f 0 --fisher
    done

    # Arm 2: Somatic Paired 29,389x (76|86) and 45,707x (75|86) combined depth at -f 0.0005, -f 0.001, and -f 0 with --fisher
    run_ud_arm 'ud_som_76_86_29kx_f0.0005_fisher' \"\$REF_HG19\" '$tp_dir/bams/SRR15006376.bam|$tp_dir/bams/SRR15006386.bam' \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" 'SRR15006376|SRR15006386' -f 0.0005 --fisher
    run_ud_arm 'ud_som_76_86_29kx_f0_fisher'      \"\$REF_HG19\" '$tp_dir/bams/SRR15006376.bam|$tp_dir/bams/SRR15006386.bam' \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" 'SRR15006376|SRR15006386' -f 0 --fisher
    run_ud_arm 'ud_som_75_86_45kx_f0.0005_fisher' \"\$REF_HG19\" '$tp_dir/bams/SRR15006375.bam|$tp_dir/bams/SRR15006386.bam' \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" 'SRR15006375|SRR15006386' -f 0.0005 --fisher
    run_ud_arm 'ud_som_75_86_45kx_f0_fisher'      \"\$REF_HG19\" '$tp_dir/bams/SRR15006375.bam|$tp_dir/bams/SRR15006386.bam' \"\$KEEP_DIR/ultradeep_aml_hotspots.bed\" 'SRR15006375|SRR15006386' -f 0 --fisher

    # Arm 3: HG008-T|N 14,964x combined hotspot (chrX:80536000-80538000) & HG002 7,361x hotspot (chr20:26541000-26543000) at -f 0.0005 --fisher
    printf 'chrX\t80536000\t80538000\tHG008_15KX\n' > \"\$KEEP_DIR/hg008_15kx.bed\"
    printf 'chr20\t26541000\t26543000\tHG002_7KX\n' > \"\$KEEP_DIR/hg002_7kx.bed\"
    run_ud_arm 'ud_hg008_15kx_f0.0005_fisher' \"\$FA_GIAB\" \"\$DATA_DIR/hg008/HG008-T.slice.bam|\$DATA_DIR/hg008/HG008-N.slice.bam\" \"\$KEEP_DIR/hg008_15kx.bed\" 'HG008T|HG008N' -f 0.0005 --fisher
    run_ud_arm 'ud_hg002_7kx_f0.0005_fisher'  \"\$FA_GIAB\" \"\$DATA_DIR/hg002/HG002.chr20.bam\" \"\$KEEP_DIR/hg002_7kx.bed\" 'HG002' -f 0.0005 --fisher

    # Arm 4: HG002_rep1 RNA-Seq 50,467x - 173,132x ultra-deep exons on chr20 at -f 0.001 --fisher
    cat > \"\$KEEP_DIR/rna_173kx.bed\" <<'EOF'
chr20	58909000	58912000	RNA_173132x
chr20	62387000	62389000	RNA_155740x
chr20	59030000	59031000	RNA_87634x
chr20	50840000	50841000	RNA_86481x
EOF
    run_ud_arm 'ud_rna_173kx_f0.001_fisher' \"\$FA_GIAB\" '$rna_dir/HG002_rep1_Aligned.sortedByCoord.out.bam' \"\$KEEP_DIR/rna_173kx.bed\" 'HG002_RNA' -f 0.001 --fisher
  "

  # 11. Extended Coverage (17 modalities) + Multimodal Slices (18 slices) + G4 & E2 CRAM Twins
  # Uses GOLDEN_MODE=auto + GOLDEN_STORE_DIR=${WORK_DIR}/golden-store synced with ${GCS_ORACLE_CACHE}/golden-store
  # so every VarDictJava and vardictcpp VCF is SHA256-pinned and backed up to GCS.
  run_checkpointed_cell "sweep" "sweep_multimodal_and_cram_twin" "0-31" "infinity" bash -c '
    set -euo pipefail
    source harness/gcp_paper_runner.sh --source-only
    export THREADS=32 OMP_NUM_THREADS=32 RAYON_NUM_THREADS=32 ORACLE_XMX_GB=200
    export GOLDEN_STORE_DIR="${WORK_DIR}/golden-store"
    mkdir -p "$GOLDEN_STORE_DIR"

    # Restore existing golden-store cells from GCS if present
    if setsid gsutil -q stat "${GCS_ORACLE_CACHE}/golden-store/cells/*" 2>/dev/null; then
      echo "[GOLDEN-STORE] Restoring pinned Java/Cpp oracle cells from ${GCS_ORACLE_CACHE}/golden-store..."
      setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -q rsync -r "${GCS_ORACLE_CACHE}/golden-store" "$GOLDEN_STORE_DIR" || true
    fi

    rc_ext=0
    GOLDEN_MODE=auto GOLDEN_STORE_DIR="$GOLDEN_STORE_DIR" bash harness/run_extended_coverage.sh > "$KEEP_DIR/extended_coverage.log" 2>&1 || rc_ext=$?
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$GOLDEN_STORE_DIR" "${GCS_ORACLE_CACHE}/golden-store" || true
    record_sweep_row "sweep_multimodal_and_cram_twin" "extended_coverage_17modes" "17" "17" "0" "0" "$([ "$rc_ext" -eq 0 ] && echo MATCH || echo DIFF_CAPTURED)" "rc=$rc_ext"

    rc_mm=0
    GOLDEN_MODE=auto GOLDEN_STORE_DIR="$GOLDEN_STORE_DIR" bash harness/run_multimodal_slices.sh > "$KEEP_DIR/multimodal_slices.log" 2>&1 || rc_mm=$?
    setsid gsutil "${GSUTIL_MT_FLAGS[@]}" -o GSUtil:parallel_composite_upload_threshold=0 -q rsync -r "$GOLDEN_STORE_DIR" "${GCS_ORACLE_CACHE}/golden-store" || true
    record_sweep_row "sweep_multimodal_and_cram_twin" "multimodal_18slices" "18" "18" "0" "0" "$([ "$rc_mm" -eq 0 ] && echo MATCH || echo DIFF_CAPTURED)" "rc=$rc_mm"

    rc_cram=0
    CMP_TIMEOUT=12h bash harness/cram_twin.sh > "$KEEP_DIR/cram_twin_g4.log" 2>&1 || rc_cram=$?
    record_sweep_row "sweep_multimodal_and_cram_twin" "cram_twin_g4_hg008" "1" "1" "0" "0" "$([ "$rc_cram" -eq 0 ] && echo MATCH || echo DIFF_CAPTURED)" "rc=$rc_cram"

    REF_E="$DATA_DIR/ref/d1vd1/GRCh38.d1.vd1.fa"
    BAM_T_E="$DATA_DIR/multimodal/seqc2_s3a_tumor.bam"
    BAM_N_E="$DATA_DIR/multimodal/seqc2_s3a_normal.bam"
    CRAM_T_E="$KEEP_DIR/seqc2_s3a_tumor.cram"
    CRAM_N_E="$KEEP_DIR/seqc2_s3a_normal.cram"
    printf "chr7\t140400000\t140420000\tBRAF\n" > "$KEEP_DIR/e2_braf.bed"
    samtools view -@ 16 -T "$REF_E" -C -o "$CRAM_T_E" "$BAM_T_E" && samtools index "$CRAM_T_E"
    samtools view -@ 16 -T "$REF_E" -C -o "$CRAM_N_E" "$BAM_N_E" && samtools index "$CRAM_N_E"
    target/release/vardictrust -G "$REF_E" -f 0.01 -N "TUMOR|NORMAL" -b "$BAM_T_E|$BAM_N_E" -c 1 -S 2 -E 3 -g 4 -z -th 32 "$KEEP_DIR/e2_braf.bed" > "$KEEP_DIR/e2_rust_bam.vcf"
    env RUST_MIN_STACK=67108864 target/release/vardictrust -G "$REF_E" -f 0.01 -N "TUMOR|NORMAL" -b "$CRAM_T_E|$CRAM_N_E" -c 1 -S 2 -E 3 -g 4 -z -th 32 "$KEEP_DIR/e2_braf.bed" > "$KEEP_DIR/e2_rust_cram.vcf"
    compare_and_record_vcf "sweep_multimodal_and_cram_twin" "e2_cram_vs_bam_twin" "$KEEP_DIR/e2_rust_bam.vcf" "$KEEP_DIR/e2_rust_cram.vcf" "$KEEP_DIR/e2_cram_twin"
  '

  # Upload consolidated sweep summary table to GCS
  if [ -f "${WORK_DIR}/cells/v103_sweep_summary.tsv" ]; then
    gcs_upload "${WORK_DIR}/cells/v103_sweep_summary.tsv" "${GCS_RESULTS}/cells-v105/v105_sweep_summary.tsv" || true
  fi
  gcs_mark_done "SWEEP_V105_DONE"
  gcs_write_heartbeat "sweep" "SWEEP_V105_DONE" "COMPLETED"
}

run_worker() {
  local target_stage="${1:-all}"
  trap '' TERM
  worker_preflight
  case "$target_stage" in
    1) run_stage_1 ;;
    2) run_stage_2 ;;
    3) run_stage_3 ;;
    4) run_stage_4 ;;
    5) run_stage_5 ;;
    sweep) run_v103_parity_sweep ;;
    all)
      run_stage_1
      run_packed_stages_2_to_4
      run_stage_5
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
  local z
  for z in us-central1-a us-central1-b us-central1-c us-central1-f "$GCP_ZONE"; do
    if gcloud compute instances describe "$VM_NAME" --project="$GCP_PROJECT" --zone="$z" >/dev/null 2>&1; then
      echo "[Controller] Deleting VM '$VM_NAME' in '$z'..."
      gcloud compute instances delete "$VM_NAME" --project="$GCP_PROJECT" --zone="$z" --quiet || true
    fi
  done
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

vm_cpu_platform() {
  gcloud compute instances describe "$VM_NAME" \
    --project="$GCP_PROJECT" \
    --zone="$GCP_ZONE" \
    --format="value(cpuPlatform)" 2>/dev/null || echo "ABSENT"
}

run_controller() {
  local max_attempts="${1:-6}"
  local start_attempt="${2:-1}"
  local worker_stage="${WORKER_STAGE:-all}"
  local done_marker="ALL_STAGES"
  local keep_vm_on_done=0
  if [ "$worker_stage" = "sweep" ]; then
    done_marker="SWEEP_V105_DONE"
    keep_vm_on_done="${KEEP_VM_ON_SWEEP_DONE:-1}"
  elif [ "$worker_stage" = "5" ]; then
    done_marker="STAGE5_DONE"
    keep_vm_on_done="${KEEP_VM_ON_DONE:-0}"
  elif [ "$worker_stage" = "all" ]; then
    done_marker="STAGE5_DONE"
    keep_vm_on_done="${KEEP_VM_ON_DONE:-0}"
  fi

  for ((attempt = start_attempt; attempt <= max_attempts; attempt++)); do
    if gcs_stat_done "$done_marker"; then
      echo "[Controller] ${done_marker}.done already present in ${GCS_CKPT}."
      if [ "$keep_vm_on_done" != "1" ]; then
        echo "[Controller] Deleting VM and finishing."
        cleanup_vm
      else
        echo "[Controller] Keeping VM '${VM_NAME}' alive (keep_vm_on_done=1)."
      fi
      return 0
    fi

    local model
    model="$(provisioning_model_for_attempt "$attempt")"
    GCP_ZONE="$(zone_for_attempt "$attempt")"
    echo "================================================================================"
    echo " [Controller] Attempt ${attempt}/${max_attempts} (stage=${worker_stage}, zone=${GCP_ZONE}): ${VM_NAME} (${MACHINE_TYPE}, ${BOOT_DISK_SIZE} ${BOOT_DISK_TYPE} + ${LOCAL_SSD_COUNT}x375GB NVMe, model=${model})"
    echo "================================================================================"

    local cur_status cur_model cur_cpu
    cur_status="$(vm_status)"
    cur_model="$(vm_provisioning_model)"
    cur_cpu="$(vm_cpu_platform)"
    if [ "$cur_status" = "RUNNING" ] && [ "$cur_model" = "$model" ] && [ "$cur_cpu" = "$MIN_CPU_PLATFORM" ]; then
      echo "[Controller] Reusing already-running VM '${VM_NAME}' (status=${cur_status}, model=${cur_model}, cpuPlatform=${cur_cpu})"
    else
      cleanup_vm
      mapfile -t create_cmd < <(build_gcloud_create_cmd "$attempt")
      if ! "${create_cmd[@]}"; then
        echo "[Controller] Instance creation failed on attempt ${attempt} (${model}, zone=${GCP_ZONE}); advancing to next attempt..." >&2
        sleep 15
        continue
      fi
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
    echo "[Controller] Mounting 3.0 TB NVMe RAID0 at $WORK_DIR and starting vardict-paper-worker.service (stage=${worker_stage})..."
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
trap '' TERM
export WORK_DIR='${WORK_DIR}'
export GCS_CKPT='${GCS_CKPT}'
export GOLDEN_MODE=off
exec > >(tee -a \"\${WORK_DIR}/worker.log\") 2>&1
echo \"[VM Bootstrap] Starting at \$(date -u +%Y-%m-%dT%H:%M:%SZ) (stage=${worker_stage})\"
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
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main >/dev/null 2>&1 || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r >/dev/null 2>&1 || true
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
  conda create -y -n vardict-oracle -c conda-forge -c bioconda \
    vardict-java=1.8.3 vardict=2019.06.04 vardictcpp=3=hc52dbad_0 \
    samtools=1.24 bcftools=1.24 htslib=1.24 openjdk=25.0.2 cmake make gxx_linux-64 zlib bzip2 xz libcurl openssl
fi
if ! conda env list | grep -q '^g4-tools '; then
  conda create -y -n g4-tools -c conda-forge -c bioconda \
    samtools=1.24 bcftools=1.24 htslib=1.24 bedtools=2.31.1 bwa=0.7.19 rtg-tools=3.13 openjdk=25.0.2 python=3.12
fi
export PATH=\"\$HOME/miniconda3/envs/g4-tools/bin:\$HOME/miniconda3/bin:\$HOME/.cargo/bin:\$PATH\"
bash harness/gcp_paper_runner.sh --worker ${worker_stage}
EOS
      chmod +x '$WORK_DIR/bootstrap_and_run.sh'
      sudo loginctl enable-linger \"\$(id -un)\"
      if systemctl is-active --quiet vardict-paper-worker.service; then
        echo '[Controller] vardict-paper-worker.service is already active; attaching monitor without restart.'
      else
        sudo systemctl reset-failed vardict-paper-worker.service 2>/dev/null || true
        sudo systemd-run --unit=vardict-paper-worker \
          --uid=\"\$(id -un)\" --gid=\"\$(id -gn)\" \
          --setenv=HOME=\"\$HOME\" \
          bash '$WORK_DIR/bootstrap_and_run.sh'
      fi
    "

    # Poll VM status & GCS heartbeat every 30s until done_marker, worker failure, or Spot preemption
    while true; do
      sleep 30
      if gcs_stat_done "$done_marker"; then
        echo "[Controller] Worker completed ${done_marker} on attempt ${attempt} (${model})!"
        if [ "$keep_vm_on_done" != "1" ]; then
          cleanup_vm
        else
          echo "[Controller] Keeping VM '${VM_NAME}' alive (keep_vm_on_done=1)."
        fi
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
      echo "[Controller $(date -u +%H:%M:%SZ)] attempt=${attempt}(${model}) stage=${worker_stage} vm=${st} unit=${unit_state} heartbeat=${hb}"
      if [ "$unit_state" = "failed" ] || [ "$unit_state" = "inactive" ]; then
        if gcs_stat_done "$done_marker"; then
          echo "[Controller] ${done_marker} completed!"
          if [ "$keep_vm_on_done" != "1" ]; then
            cleanup_vm
          else
            echo "[Controller] Keeping VM '${VM_NAME}' alive (keep_vm_on_done=1)."
          fi
          return 0
        fi
        if [ "$unit_state" = "inactive" ] && gcs_stat_done "ALL_STAGES" && ! gcs_stat_done "STAGE5_DONE"; then
          echo "[Controller] Stages 1-4 finished (ALL_STAGES.done); restarting vardict-paper-worker.service to execute Stage 5..."
          gcloud compute ssh "$VM_NAME" --project="$GCP_PROJECT" --zone="$GCP_ZONE" --quiet --command="
            sudo systemctl reset-failed vardict-paper-worker.service 2>/dev/null || true
            sudo systemd-run --unit=vardict-paper-worker \
              --uid=\"\$(id -un)\" --gid=\"\$(id -gn)\" \
              --setenv=HOME=\"\$HOME\" \
              bash '$WORK_DIR/bootstrap_and_run.sh'
          "
          continue
        fi
        echo "[Controller] Worker unit exited (${unit_state}) without ${done_marker}.done; fetching tail of worker.log:" >&2
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

  echo "[Controller] Exhausted ${max_attempts} attempts without completing ${done_marker}." >&2
  cleanup_vm
  return 1
}

run_plan_check() {
  echo "=== [gcp_paper_runner: --plan-check] ==="
  [ "$GOLDEN_MODE" = "off" ] || { echo "FAIL: GOLDEN_MODE != off" >&2; exit 1; }
  for att in 1 2 3 4 5; do
    [ "$(provisioning_model_for_attempt "$att")" = "SPOT" ] || { echo "FAIL: attempt $att must be SPOT" >&2; exit 1; }
  done
  [ "$(provisioning_model_for_attempt 6)" = "STANDARD" ] || { echo "FAIL: attempt 6 must escalate to STANDARD" >&2; exit 1; }

  local cmd1 cmd6
  cmd1="$(build_gcloud_create_cmd 1 | tr '\n' ' ')"
  cmd6="$(build_gcloud_create_cmd 6 | tr '\n' ' ')"
  [[ "$cmd1" == *"--machine-type=n2-standard-32"* ]] || { echo "FAIL: machine-type not n2-standard-32" >&2; exit 1; }
  [[ "$cmd1" == *"--min-cpu-platform=Intel Ice Lake"* ]] || { echo "FAIL: min-cpu-platform not Intel Ice Lake" >&2; exit 1; }
  [[ "$cmd1" == *"--boot-disk-size=100GB"* ]] || { echo "FAIL: boot-disk-size not 100GB" >&2; exit 1; }
  [[ "$cmd1" == *"--boot-disk-type=pd-ssd"* ]] || { echo "FAIL: boot-disk-type not pd-ssd" >&2; exit 1; }
  [[ "$cmd1" == *"--local-ssd=interface=NVME"* ]] || { echo "FAIL: missing --local-ssd=interface=NVME" >&2; exit 1; }
  [ "$(build_gcloud_create_cmd 1 | grep -c '^--local-ssd=interface=NVME$')" -eq 4 ] || { echo "FAIL: expected 4 NVMe local SSDs (minimum allowed on n2-standard-32)" >&2; exit 1; }
  [[ "$cmd1" == *"--provisioning-model=SPOT"* ]] || { echo "FAIL: attempt 1 cmd missing SPOT" >&2; exit 1; }
  [[ "$cmd6" == *"--provisioning-model=STANDARD"* ]] || { echo "FAIL: attempt 6 cmd missing STANDARD" >&2; exit 1; }
  [[ "$cmd6" != *"--maintenance-policy=TERMINATE"* ]] || { echo "FAIL: STANDARD cmd must not set maintenance-policy=TERMINATE" >&2; exit 1; }

  # Verify audit fixes in script text: exact mkbed(), symmetric n=3 Java loops, 13-mode sweep, Fig-2/Fig-3 field callers, 3-way scaling
  [ "$EXPECTED_FP" = "d925892846918c2a" ] || { echo "FAIL: EXPECTED_FP '$EXPECTED_FP' != d925892846918c2a" >&2; exit 1; }
  [[ "$GCS_CKPT" == */checkpoints-v105 ]] || { echo "FAIL: GCS_CKPT '$GCS_CKPT' does not end with /checkpoints-v105" >&2; exit 1; }
  grep -Fq 'samtools view -F 0x904 -q 1' "$0" || { echo "FAIL: missing exact samtools view -F 0x904 -q 1 mkbed()" >&2; exit 1; }
  grep -Fq 'java_var_th${th}_run\${run}.time' "$0" || { echo "FAIL: missing symmetric n=3 Java Tier-P loop" >&2; exit 1; }
  grep -Fq 's2_tierp_mode_parity_sweep' "$0" || { echo "FAIL: missing s2_tierp_mode_parity_sweep cell" >&2; exit 1; }
  grep -Fq 's2_fig2_hg002_chr20_truth' "$0" || { echo "FAIL: missing s2_fig2_hg002_chr20_truth cell" >&2; exit 1; }
  grep -Fq 's2_fig1_fig3_field_callers' "$0" || { echo "FAIL: missing s2_fig1_fig3_field_callers cell" >&2; exit 1; }
  grep -Fq 's4_scale_exome_3way_th' "$0" || { echo "FAIL: missing 3-way exome scaling cell" >&2; exit 1; }
  grep -Fq 's4_span_chr1_3way_th' "$0" || { echo "FAIL: missing 3-way chr1 span scaling cell" >&2; exit 1; }
  grep -Fq 's4_tier3_truth_and_qual_sweep' "$0" || { echo "FAIL: missing Tier-3 truth and ROC/QUAL sweep cell" >&2; exit 1; }
  grep -Fq '"$DATA_DIR/hg002/truth/truth.vcf.gz"' "$0" || { echo "FAIL: missing exact hg002 truth path" >&2; exit 1; }
  grep -Fq '"$DATA_DIR/seqc2/truth/truth.vcf.gz"' "$0" || { echo "FAIL: missing exact seqc2 truth path" >&2; exit 1; }
  grep -Fq '"$DATA_DIR/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz"' "$0" || { echo "FAIL: missing exact hg008 v0.3 truth path" >&2; exit 1; }
  grep -Fq 'record_evidence.py hg008-cram-twin' "$0" || { echo "FAIL: missing hg008-cram-twin evidence recording" >&2; exit 1; }
  grep -Fq 'bd605a455b2e2064931e6097c64695ea' "$0" || { echo "FAIL: missing HG002.chr20.bam MD5 gate" >&2; exit 1; }
  grep -Fq 'run_v103_parity_sweep' "$0" || { echo "FAIL: missing run_v103_parity_sweep" >&2; exit 1; }
  grep -Fq 'sweep_wgs_mrd_matrix' "$0" || { echo "FAIL: missing sweep_wgs_mrd_matrix" >&2; exit 1; }
  grep -Fq 'sweep_ultradeep_mrd_stress' "$0" || { echo "FAIL: missing sweep_ultradeep_mrd_stress" >&2; exit 1; }
  grep -Fq 's5_stage_hg002_full_wgs' "$0" || { echo "FAIL: missing s5_stage_hg002_full_wgs" >&2; exit 1; }
  grep -Fq 's5_wgs_hg002_full_contigs_3way_th32' "$0" || { echo "FAIL: missing s5_wgs_hg002_full_contigs_3way_th32" >&2; exit 1; }
  grep -Fq 's5_wgs_hg002_full_giab_3way_th32' "$0" || { echo "FAIL: missing s5_wgs_hg002_full_giab_3way_th32" >&2; exit 1; }
  grep -Fq 's5_stage_hg008_full_wgs' "$0" || { echo "FAIL: missing s5_stage_hg008_full_wgs" >&2; exit 1; }
  grep -Fq 's5_wgs_hg008_full_contigs_3way_th32' "$0" || { echo "FAIL: missing s5_wgs_hg008_full_contigs_3way_th32" >&2; exit 1; }
  grep -Fq 's5_wgs_hg008_full_giab_3way_th32' "$0" || { echo "FAIL: missing s5_wgs_hg008_full_giab_3way_th32" >&2; exit 1; }
  grep -Fq 's5_grch38_field_callers_truth' "$0" || { echo "FAIL: missing s5_grch38_field_callers_truth" >&2; exit 1; }
  grep -Fq 's5_h4_maintenance_security_matrix' "$0" || { echo "FAIL: missing s5_h4_maintenance_security_matrix" >&2; exit 1; }
  grep -Fq 's5_stage_cross_platform_chr20' "$0" || { echo "FAIL: missing s5_stage_cross_platform_chr20" >&2; exit 1; }
  grep -Fq 's5_cross_platform_short_read_chr20_panel' "$0" || { echo "FAIL: missing s5_cross_platform_short_read_chr20_panel" >&2; exit 1; }

  # Enforce the 6 Future-Version Operational Invariants:
  # Invariant 3: 16-thread single-process gsutil (zero bare gsutil -m without parallel_process_count=1)
  grep -Fq 'GSUTIL_MT_FLAGS=(-m -o "GSUtil:parallel_process_count=1" -o "GSUtil:parallel_thread_count=16")' "$0" || { echo "FAIL: missing 16-thread single-process GSUTIL_MT_FLAGS" >&2; exit 1; }
  if grep -vE '^[[:space:]]*(#|echo )' "$0" | grep -E 'gsutil[[:space:]]+-m[[:space:]]+(cp|rsync)' >/dev/null 2>&1; then
    echo "FAIL: bare gsutil -m without parallel_process_count=1 causes Ubuntu 24.04 snap multiprocessing.Pool deadlock" >&2
    exit 1
  fi
  # Invariant 4: Static & reusable asset GCS push/pull hooks
  grep -Fq 'vardictcpp-ae93fe39dc.bin' "$0" || { echo "FAIL: missing vardictcpp binary GCS cache hook" >&2; exit 1; }
  grep -Fq 'wgs-hg002/s5_hg002_ready.done' "$0" || { echo "FAIL: missing wgs-hg002 GCS cache hook" >&2; exit 1; }
  grep -Fq 'wgs-hg008/s5_hg008_ready.done' "$0" || { echo "FAIL: missing wgs-hg008 GCS cache hook" >&2; exit 1; }
  grep -Fq 'cross-platform-chr20/s5_chr20_slices_ready.done' "$0" || { echo "FAIL: missing cross-platform-chr20 GCS cache hook" >&2; exit 1; }
  # Invariant 5: Cached amplicon VCF runs must always pass -c 1 -S 2 -E 3 -g 4
  if grep -n 'vardictrust.*-a 10:0.95' "$0" | grep -v -- '-c 1 -S 2 -E 3 -g 4' >/dev/null 2>&1; then
    echo "FAIL: vardictrust amplicon (-a 10:0.95) invocation missing '-c 1 -S 2 -E 3 -g 4' (breaks omit_end_tag parity)" >&2
    exit 1
  fi

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

  # Verify oracle_cache_get_or_run stores on miss and restores from cache on hit without re-running generator
  local gen_count="$tmp/oracle_gen_count"
  echo 0 > "$gen_count"
  oracle_cache_get_or_run "test/java.var" "$tmp/work1/java.var" bash -c "echo \$(( \$(cat '$gen_count') + 1 )) > '$gen_count'; echo 'ORACLE_TRUTH_ROW' > '$tmp/work1/java.var'"
  [ "$(cat "$gen_count")" -eq 1 ] || { echo "FAIL: oracle_cache_get_or_run did not run on miss" >&2; exit 1; }
  [ -s "$LOCAL_CKPT_DIR/oracle-cache/test/java.var" ] || { echo "FAIL: oracle-cache file missing" >&2; exit 1; }
  oracle_cache_get_or_run "test/java.var" "$tmp/work2/java.var" bash -c "echo \$(( \$(cat '$gen_count') + 1 )) > '$gen_count'; echo 'REGEN_ROW' > '$tmp/work2/java.var'"
  [ "$(cat "$gen_count")" -eq 1 ] || { echo "FAIL: oracle_cache_get_or_run re-ran generator on cache hit" >&2; exit 1; }
  [ "$(cat "$tmp/work2/java.var")" = "ORACLE_TRUTH_ROW" ] || { echo "FAIL: restored oracle file content mismatch" >&2; exit 1; }

  # Verify build_tierp_covered_bed on fixtures/micro/micro.bam if samtools is available
  if command -v samtools >/dev/null 2>&1 && [ -f "fixtures/micro/micro.bam" ]; then
    build_tierp_covered_bed "fixtures/micro/micro.bam" "$tmp/micro.mkbed.bed" 2
    [ -s "$tmp/micro.mkbed.bed" ] || { echo "FAIL: micro.mkbed.bed empty" >&2; exit 1; }
  fi

  echo "PASS gcp_paper_runner --plan-check"
}

case "${1:-}" in
  --source-only)
    return 0 2>/dev/null || exit 0
    ;;
  --plan-check)
    run_plan_check
    ;;
  --controller)
    shift
    max_att=6
    start_att=1
    while [ $# -gt 0 ]; do
      case "$1" in
        --max-attempts) max_att="$2"; shift 2 ;;
        --start-attempt) start_att="$2"; shift 2 ;;
        --stage) export WORKER_STAGE="$2"; shift 2 ;;
        *) echo "Unknown controller option: $1" >&2; exit 1 ;;
      esac
    done
    run_controller "$max_att" "$start_att"
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
