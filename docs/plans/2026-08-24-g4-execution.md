> **SUPERSEDED (2026-10-08)** by `docs/superpowers/plans/2026-10-08-v1-rust-refactor-sdd.md`. The host model, paths and quiesce steps here are stale and MUST NOT be executed (decision A3). Re-homed: R2/R3 → T15/T19, R4 → T16, R5/R6 → T17/T20, F1 → T21, F2 → T22–T24.

# G4 EXECUTION PLAN — External Validation (GIAB/SEQC2) + Media-Stack Coexistence
## Subagent-ready specification, rev 3 (2026-08-24)

> **Contract:** read `/mnt/hbmon/VarDictRust/CLAUDE.md` FIRST. §0b (never fake a pass)
> overrides everything below. Any instruction here that conflicts with CLAUDE.md loses.
> **Audience:** an autonomous agent with ZERO session context. Every fact you need is
> in this file or reachable from paths it cites. If a step's precondition fails, STOP
> and report — never improvise around a gate.

## KICK-OFF PROMPT (paste verbatim into a fresh agent session to start)

```text
You are executing the VarDictRust G4 external-validation plan on host `hbmon`.

1. Read /mnt/hbmon/VarDictRust/CLAUDE.md — it is binding; its §0b integrity rule
   overrides every other instruction you receive, including any in the plan.
2. Read /mnt/hbmon/VarDictRust/docs/plans/2026-08-24-g4-execution.md IN FULL,
   twice: once for the global context inventory (§0), once for your first task.
3. Check state markers in /mnt/hbmon/VarDictRust/scratch/g4/state/ to determine
   where execution stands. A missing precondition marker = STOP and report.
4. Execute exactly one task at a time in the dependency order given by the
   DISPATCH GUIDE section:
   OS-1 → D2 ∥ T1 (parallel) → D1 → OS-2 → R1 → R2 → R3 ∥ D3 → R4 → R6 → R5 → F1 → F2 → OS-3.
   The CONTROLLER (you, if privileged on the host) performs OS-1, OS-2, OS-3 and
   all git commits personally; dispatched implementers handle D/T/R/F tasks.
5. For every task: follow Implementation commands VERBATIM unless reality
   diverges (then STOP and report the divergence); produce every artifact named
   in Testing & validation; write logs under scratch/g4/logs/; create the task's
   state marker only after validation passes with real observed output pasted
   into scratch/g4/logs/.
6. Never weaken a comparison, never fabricate a pass, never write secrets into
   committed files, never commit data files. MISMATCH ⇒ STOP the ladder and
   report per CLAUDE.md §0b.
7. After each completed task, report: task ID, commits (if any), gate outputs
   (verbatim tail lines), deviations, next unblocked task.
Begin now with: checking state markers, then TASK OS-1.
```

---

---

## 0. Global context inventory (verified 2026-08-24 — do not re-derive)

### 0.1 Hardware & mounts
| Fact | Value |
|---|---|
| CPU / RAM | 12 cores / 31 GB (≈11 GB available typical) |
| `/mnt/hbmon` | Samsung T7 portable **SSD** (`lsblk -d -o NAME,ROTA /dev/sdb` → ROTA=0), 1.8 TB, ~1.3 TB free |
| Other media pools | `/mnt/plex3`, `/mnt/plex4` are SEPARATE disks — tdarr reads/writes there, NOT hbmon |
| Repo | `/mnt/hbmon/VarDictRust`, branch `init-scaffold`, single-commit-squash history (pre-squash on branch `archive/pre-squash-20260823`) |
| Scratch root | `/mnt/hbmon/VarDictRust/scratch/g4/` (repo `.gitignore` covers `scratch/`; create freely) |

### 0.2 VarDict toolchain (already working — do not rebuild)
- Our binary: build via `cargo build --release` → `target/release/vardictrust`
- Oracle: conda env **`vardict-oracle`** wrapping VarDictJava v1.8.3 (@`0886261`)
- Harness (read each script's header BEFORE first use):
  - `harness/run_oracle.sh <ref.fa> <bam-or-pair> <bed> <out.vcf> [sample]` — runs Java engine; somatic pair syntax `"tumor.bam|normal.bam"`
  - `harness/compare.sh <ref.fa> <bam-or-pair> <bed> [sample]` — dual-engine Tier-2 VCF compare; **exit 0 + prints MATCH = pass**
  - `harness/slice.sh <in.bam> <chr:start-end> <out.bam>` — slice big BAMs
  - `bash harness/run_integration.sh` — full corpus; healthy state = **91 PASS / 1 FAIL** (fail = `hard_clip`, permanent stale fixture). Run once before starting (record output) and after ANY change touching engine/harness.

### 0.3 Docker media stack (compose project: `/home/palak/Documents/docker-compose.yml`)
Containers relevant here: `deluge` (:64 in compose), `tdarr` (:267), `sabnzbd` (:537, host port **8081**→8080). All have `restart=always` (irrelevant to manual stops). `vpn` container serves deluge's chain — **never touch it**.

### 0.4 Host automation that WILL fight you (all code-verified 2026-08-24)
| Controller | Mechanism (file:line) | Threat | Handling |
|---|---|---|---|
| `~/.config/systemd/user/media-stream-qos.service` → `media_stream_qos_daemon.py` | Auto-**resumes SAB** after 45 s idle-debounce (`:273-275`); throttles/unthrottles Deluge (`:247`); pauses tdarr nodes under CPU load (`:286`); writes heartbeat `/dev/shm/media_streaming_state.json` | Reverts our queue pauses ≤1 loop | **Mask** at M0: `systemctl --user mask --now media-stream-qos.service`; unmask at M2 |
| cron `deluge_stalled_watchdog.py` (*/2 min) | resumes torrents (`:650-750`), `force_resume_seeding` (`:750`), SAB `mode=resume` (`:988`) | Undoes pauses ≤2 min | Gate cron line M0→M2 |
| cron `sab_queue_tier_engine.py` (*/5 min) | promotes (=resumes) paused slots (`:334`) despite docstring | Partial unpause ≤5 min | Gate M0→M2 |
| cron `deluge_fake_download_watchdog.py` (*/5 + daily) | audits/manages deluge state | Interference | Gate M0→M2 |
| cron `tdarr_guardian.py` (*/20 min) | bounded `docker restart tdarr` on stagnation (`:197,249-257`) | Revives stopped tdarr | Gate M0→**M3** |
| cron `tdarr_swap_wrapper.py` (*/10 min) | swaps tdarr caches, holds `/tmp/media_ops.lock` | Interference/log spam | Gate M0→M3 |
| cron `auto_torrent_picker_engine.py` (6 h) + `sonarr_season_upgrader.py --upgrade-next` (30 min) + `radarr_movie_upgrader.py --upgrade-next` (10 min) + `sonarr_missing_episodes_filler.py` (20 min) | push new download jobs into clients | Job pileup while clients down | Gate M0→M2 |
| cron `stack_nightly_maintenance.py` (03:45 daily) | per-service `docker compose stop` → SQLite VACUUM → **`docker compose start`** (`:301,329`) | Would RESURRECT stopped containers overnight | Gate entire window M0→M3 |
| cron `hbmon_disk_breaker.py` (*/5, HB_FLOOR=150 HB_TARGET=300) | purges torrents only when hbmon free <150 GB; RPC-failure-safe by design | NONE at ~1300 GB free | Leave running |
| Torznab bridges / homepage-proxy / audiobook+ABS+bazarr jobs / DB backups / price watchdog | do not touch our three containers | none | Leave running |

### 0.5 Secrets policy (BINDING)
API keys/passwords exist in `/mnt/hbmon/app_configs/*/` and inside user scripts.
**Never** write a secret value into any committed file, log, or the MANIFEST.
Read at runtime when needed, e.g.:
```bash
SAB_KEY=$(awk '$1=="api_key"{print $3}' /mnt/hbmon/app_configs/sabnzbd/sabnzbd.ini)
```

### 0.6 State protocol for parallel/dispatched tasks
All tasks coordinate via `scratch/g4/state/`. Before doing work, check markers;
after verified completion, create yours:
```
scratch/g4/state/{m0.done,m1.done,m2.done,m3.done,
                  d1-hcc1395-bams.done,d2-truth-ref.done,d3-hg008.done,t1-tools.done,
                  r1-regions.done,r2-slicegate.done,r3-fullwes.done,r4-hg008.done,
                  r5-truth.done,r6-cram.done}
```
Marker content: one-line JSON `{"task":"…","ts":"ISO8601","evidence":"path"}`.
A task whose precondition marker is missing STOPS and reports.

### 0.7 Time & disk budgets
Downloads are bandwidth-bound (est. 20–40 GB total); full-WES dual-engine ≈ hours on
12 c. Whole G4 realistically straddles overnight. Keep ≥100 GB free on hbmon at all
times (check before each large write; abort downloads if approaching).

### 0.8 Global abort criteria (any hit ⇒ STOP, report, do not improvise)
1. A `compare.sh` result is neither MATCH nor an understood documented exception → treat as MISMATCH per CLAUDE.md §0b; minimize + report. NEVER weaken comparators/eps/cols.
2. hbmon free space < 100 GB.
3. Oracle env missing/broken (`conda run -n vardict-oracle vardict-java -h` fails) → run `bash harness/setup_oracle.sh` first; report if that fails.
4. Any secret would land in a committed file.
5. Corpus baseline changed: re-running `run_integration.sh` gives anything other than 91 PASS / 1 FAIL(hard_clip) → something regressed; STOP.

---

## ORCHESTRATION TASKS (host-level; run these yourself, not in containers)

### TASK OS-1 — Media-stack quiesce (Milestone M0 + M1)
**Problem statement.** G4 needs predictable disk/CPU on hbmon and stable queue states
for hours, but six cron jobs + one systemd daemon actively pause/resume/restart the
three containers. We must quiesce them *gracefully* (let actives finish), *reversibly*
(one-command restore), and *verifiably*.

**Solution theory.** Four layers of control, applied in dependency order:
(1) mask the autonomous daemon (it owns queue state), (2) gate the cron controllers
(backup → filtered crontab → restore-from-backup later), (3) pause queues so actives
drain instead of dying mid-write, (4) `docker stop` with generous SIGTERM grace.
Reversibility comes from artifacts: a timestamped crontab backup + systemd unmask.

**Implementation.**
```bash
set -euo pipefail
cd /mnt/hbmon/VarDictRust
mkdir -p scratch/g4/state scratch/g4/logs
STAMP=$(date +%Y%m%d-%H%M%S)

# --- M0a: mask the QoS daemon (it owns pause/resume of SAB + deluge speeds)
systemctl --user mask --now media-stream-qos.service
systemctl --user is-active media-stream-qos.service || true   # expect "inactive"/"masked"

# --- M0b: gate cron controllers (keep backup; restore verbatim at M2/M3)
crontab -l > ~/Documents/crontab.backup.$STAMP
echo "$STAMP" > scratch/g4/state/crontab-backup-stamp
GATE='tdarr_guardian\.py|tdarr_swap_wrapper\.py|sab_queue_tier_engine\.py|deluge_stalled_watchdog\.py|deluge_fake_download_watchdog\.py|auto_torrent_picker_engine\.py|stack_nightly_maintenance\.py|sonarr_season_upgrader\.py --upgrade-next|radarr_movie_upgrader\.py --upgrade-next|sonarr_missing_episodes_filler\.py'
crontab -l | grep -Ev "$GATE" | crontab -
crontab -l > scratch/g4/logs/crontab.gated.$STAMP          # evidence copy

# --- M1a: graceful queue drain (actives finish; nothing new starts)
SAB_KEY=$(awk '$1=="api_key"{print $3}' /mnt/hbmon/app_configs/sabnzbd/sabnzbd.ini)
curl -sf "http://127.0.0.1:8081/api?apikey=$SAB_KEY&mode=pause"
docker exec deluge deluge-console "pause *" || echo "WARN: deluge-console pause rc=$?"

# --- M1b: wait for queues to go idle (cap 30 min; SIGTERM checkpoint is still clean)
for i in $(seq 1 60); do
  SAB_ACTIVE=$(curl -sf "http://127.0.0.1:8081/api?apikey=$SAB_KEY&mode=qstatus&output=json" \
               | python3 -c 'import json,sys;d=json.load(sys.stdin);print(sum(1 for s in d.get("slots",{}) if s.get("status") in ("Downloading","Verifying","Repairing","Extracting")))' 2>/dev/null || echo 9)
  DL_STATE=$(docker exec deluge deluge-console "info" 2>/dev/null | grep -c "Downloading" || echo 9)
  echo "poll$i sab_active=$SAB_ACTIVE deluge_downloading=$DL_STATE" >> scratch/g4/logs/drain.log
  [ "$SAB_ACTIVE" = "0" ] && [ "$DL_STATE" = "0" ] && break
  sleep 30
done

# --- M1c: spin down (order: heaviest I/O last-in first-out; generous grace)
docker stop -t 60 sabnzbd deluge
docker stop -t 180 tdarr
docker ps --format '{{.Names}}' | grep -xE 'deluge|sabnzbd|tdarr' && { echo "FAIL: still up"; exit 1; } || true
touch scratch/g4/state/m0.done scratch/g4/state/m1.done
```

**Testing & validation.**
- `crontab -l | grep -Ec 'tdarr_guardian|sab_queue_tier'` → **0**
- `systemctl --user is-enabled media-stream-qos.service` → contains `masked`
- `docker ps` shows none of the three; `vpn` still Up
- Drain log shows a trailing `sab_active=0 deluge_downloading=0` OR 30-min cap reached (either acceptable; note which in report)
- Baseline corpus sanity (once, before downloads saturate disk): `bash harness/run_integration.sh 2>&1 | tail -5` → PASS 91 / FAIL 1 hard_clip. Record in `scratch/g4/logs/corpus-baseline.txt`.

**Rollback (if aborted mid-task).** Restore crons: `crontab ~/Documents/crontab.backup.$STAMP`;
unmask daemon: `systemctl --user unmask media-stream-qos.service && systemctl --user start media-stream-qos.service`;
restart containers: `docker start tdarr sabnzbd deluge`.

---

### TASK OS-2 — Milestone M2 transition (run ONLY after D1+D2+D3 all show `.done`)
**Problem statement.** Once every vardict byte is downloaded and md5-verified, the
media stack's downloader half should return to normal duty while G4 continues
compute-only work — without letting controllers loose before we deliberately set
queue state.

**Solution theory.** Order is load-bearing: containers healthy → deliberate queue
resume → downloader crons restored → QoS daemon unmasked LAST (so it inherits our
intended state rather than imposing its own). tdarr-side gates stay on until M3.

**Implementation.**
```bash
set -euo pipefail
cd /mnt/hbmon/VarDictRust
for m in d1-hcc1395-bams d2-truth-ref d3-hg008; do
  [ -f scratch/g4/state/$m.done ] || { echo "M2 blocked: $m missing"; exit 1; }
done
STAMP=$(cat scratch/g4/state/crontab-backup-stamp)
BACKUP=~/Documents/crontab.backup.$STAMP

docker start sabnzbd deluge
sleep 20   # daemons bind ports
SAB_KEY=$(awk '$1=="api_key"{print $3}' /mnt/hbmon/app_configs/sabnzbd/sabnzbd.ini)
curl -sf "http://127.0.0.1:8081/api?apikey=$SAB_KEY&mode=resume"
docker exec deluge deluge-console "resume *" || true

# restore ONLY downloader-side cron lines (tdarr pair + nightly stay gated)
UNGATE_DL='sab_queue_tier_engine|deluge_stalled_watchdog|deluge_fake_download_watchdog|auto_torrent_picker_engine|--upgrade-next|missing_episodes_filler'
{ grep -vE 'tdarr_guardian\.py|tdarr_swap_wrapper\.py|stack_nightly_maintenance\.py' "$BACKUP"; } > /tmp/cron.m2
# re-add gated-away downloader lines from backup:
grep -E "$UNGATE_DL" "$BACKUP" >> /tmp/cron.m2
sort -u /tmp/cron.m2 | crontab -

systemctl --user unmask media-stream-qos.service
systemctl --user start  media-stream-qos.service
touch scratch/g4/state/m2.done
```

**Testing & validation.**
- `curl -sf "...mode=qstatus&output=json"` parses; SAB reports queue (paused count may be whatever tiering decides — that's normal post-restore)
- `docker ps`: deluge+sabnzbd `(healthy)` within ~2 min; tdarr still absent
- `crontab -l | grep -c tdarr_guardian` → **0** (still gated); `grep -c sab_queue_tier_engine` → **1**
- `systemctl --user is-active media-stream-qos.service` → `active`

**Rollback.** Re-mask daemon + `crontab $(ls -t ~/Documents/crontab.backup.* | head -1)` + stop the two containers.

---

### TASK OS-3 — Milestone M3 handback (after R3–R6 all `.done`)
**Problem statement.** Return the box to its owner exactly as found: tdarr running,
every cron line present, daemon live, no leftover state.

**Implementation.**
```bash
cd /mnt/hbmon/VarDictRust
STAMP=$(cat scratch/g4/state/crontab-backup-stamp)
docker start tdarr
sleep 45; docker inspect -f '{{.State.Health.Status}}' tdarr    # want "healthy"
crontab ~/Documents/crontab.backup.$STAMP                        # verbatim restore
diff <(crontab -l) ~/Documents/crontab.backup.$STAMP && echo CRON-RESTORED
touch scratch/g4/state/m3.done
```

**Testing & validation.** `docker ps` → all three `(healthy)`; diff above prints
CRON-RESTORED; `systemctl --user is-active media-stream-qos.service` → active.

---

## DATA TASKS

### TASK D1 — fetch_data.sh + SEQC2 HCC1395 WES tumor/normal BAMs
**Problem statement.** PLAN.md mandates real-data validation on SEQC2 HCC1395 WES
(T167). We need the paired tumor/normal exome BAMs (plus BAI indices) from NIST/GIAB
hosting, integrity-verified, with full provenance, reproducible by script.

**Solution theory.** NCBI GIAB hosts somatic benchmark resources under
`https://ftp.ncbi.nlm.nih.gov/giab/` (HTTPS works — verified HTTP 200). Exact asset
names drift over years, so the script DISCOVERS candidates from index pages, applies
selection rules, verifies published MD5s, and records everything in
`benchmarks/MANIFEST` (committed; contains URLs+md5s+dates+commands, never secrets).
Use HTTPS + `wget -c` (resumable); two parallel streams max (SSD + pipe friendly).

**Implementation.** Write `harness/fetch_data.sh` implementing:
1. `BASE=https://ftp.ncbi.nlm.nih.gov/giab`; browse
   `$BASE/data_somatic/` index (curl) to locate the HCC1395 directory (PLAN.md
   research cites this tree name; if renamed, follow index links and RECORD the redirect chain).
2. Select files matching: tumor BAM ≥1 GB containing `tumor|HCC1395_T|T-`; normal
   BAM containing `normal|HCC1395_N|N-`; prefer WES/exome in filename over WGS;
   take `.bai`/`.bam.bai` siblings; grab sibling `.md5` files when present.
3. Download to `scratch/g4/hcc1395/` with `wget -c -q --show-progress`; two files concurrently.
4. Verify: `md5sum -c` against fetched/published md5 (when the site provides none,
   record OUR computed md5 in MANIFEST marked `computed-local` — honest provenance).
5. Append a MANIFEST block per file: URL, final path, bytes, md5, md5-source, date,
   command used. MANIFEST lives at `benchmarks/MANIFEST`.
6. On success: `touch scratch/g4/state/d1-hcc1395-bams.done`.

**Testing & validation.**
- `md5sum -c` exit 0 for every artifact (or `computed-local` recorded honestly)
- `samtools quickcheck -v *.bam` silent (valid EOF/BAI pairing)
- Both BAMs open: `samtools view -H tumor.bam | head -3` prints an @HD/@SQ header
- MANIFEST block present per file; no secret strings anywhere (`grep -R "$SAB_KEY\|DELUGE_PASSWORD" benchmarks/ || echo clean`)
- Script is idempotent: second run skips existing verified files (`[ -f X ] && md5 ok && skip`)

**Failure handling.** URL 404 → walk parent index, record new path in MANIFEST; md5
mismatch → delete + redownload once; second mismatch → STOP & report (corrupt source).

---

### TASK D2 — Truth assets + GRCh38 reference
**Problem statement.** Region selection (R1) and truth-tier scoring (R5) need the
SEQC2 HCC1395 truth VCF + high-confidence BED, GIAB HG008 truth VCF + BED, and a
GRCh38 reference FASTA(+.fai). These are small except the reference (~3 GB) — fetch
truth FIRST so R1 can proceed while big BAMs stream.

**Solution theory.** Same discover-verify-record pattern as D1, pointed at the truth
subtrees (`data_somatic/HCC1395*`, `data_somatic/HG008/Liss_lab/` per PLAN.md §622)
and NCBI/GIAB reference sets. Reference choice: GRCh38 **without alt/decoy** analysis
set unless the truth VCF headers dictate otherwise; MUST match BAM contig naming
(`chr`-prefixed). After download: `samtools faidx` + verify a known sequence.

**Implementation.** Extend `fetch_data.sh` with `--truth` and `--ref` modes writing to
`scratch/g4/truth/` and `scratch/g4/ref/`. Steps mirror D1 (discover → select
highest-version VCF/BED pairs → wget → md5 → MANIFEST). Reference: locate GRCh38
FASTA on GIAB/NCBI (e.g., `reference_materials/` subtrees), download + `samtools faidx`,
spot-check: extract chr1:1-60 twice → identical output.
Success markers: `d2-truth-ref.done`.

**Testing & validation.**
- `bcftools view -h truth.vcf | grep -i "^##contig" | head -2` shows `chr`-prefixed contigs consistent with BAM headers (`samtools view -H tumor.bam | grep ^@SQ | head -3`)
- `bcftools index -t truth.vcf.gz` succeeds (index it if gzipped)
- BED intervals sorted, within contig lengths (`bedtools`-free check: awk monotonic start per chrom)
- FAI exists; `samtools faidx ref.fa chr1:1-60` twice → byte-equal
- MANIFEST blocks complete

---

### TASK D3 — GIAB HG008 somatic BAM(s) (T168 inputs)
**Problem statement.** Second dataset (PDAC somatic) for T168 slice validation.
**Solution theory.** Identical to D1 against the HG008/Liss_lab subtree; slices only
are used downstream, but fetch the full source BAM once (slicing is local).
**Implementation.** As D1, into `scratch/g4/hg008/`; marker `d3-hg008.done`.
**Testing & validation.** md5 + `samtools quickcheck` + MANIFEST, as D1. May run
CONCURRENTLY with R3 compute (I/O-light vs compute-heavy; SSD tolerates both).

---

## TOOLING TASK

### TASK T1 — Analysis tools env
**Problem statement.** R5 needs `rtg-tools` (vcfeval) and `truvari`; R1/R4 need
`bcftools`. Isolate from the oracle env.
**Solution theory.** One dedicated conda env keeps versions pinned in MANIFEST and
keeps `vardict-oracle` pristine.
**Implementation.**
```bash
conda create -y -n g4-tools -c bioconda -c conda-forge rtg-tools truvari bcftools samtools
conda run -n g4-tools rtg version && conda run -n g4-tools truvari --version
conda run -n g4-tools bcftools --version | head -1
```
Record versions into MANIFEST. Marker `t1-tools.done`. Can run during M0/M1 drain.
**Validation.** All three version prints succeed; `rtg format` help accessible.
**Rollback.** `conda env remove -n g4-tools` (nothing else depends on it).

---

## ANALYSIS TASKS (the actual science)

### TASK R1 — Slice-region selection
**Problem statement.** Validate cheaply before committing to full-WES: pick 3 exome
regions per dataset that contain truth-set variants (so a working caller produces
non-empty output), sized minutes-per-run.
**Solution theory.** Intersect truth VCF positions with the high-confidence BED,
cluster nearby variants into windows (±2 kb padding, merge overlaps), rank by variant
count, take top-3 well-separated (>5 Mb apart) windows. Deterministic, documented.
**Implementation.**
```bash
source activate g4-tools 2>/dev/null || conda activate g4-tools
bcftools view -R highconf.bed -V indels truth.vcf.gz -Ou | bcftools query -f '%CHROM\t%POS\n' > scratch/g4/r1/snvs.txt
# cluster: sort, pad ±2000, merge overlaps, drop windows <5 variants, take top-3 spread ≥5Mb
# emit scratch/g4/r1/regions.bed  (tab: CHROM START END NAME) + human-readable rationale in r1/README.md
```
(Exact awk/sort pipeline goes in `harness/fetch_data.sh --make-regions` or a small
committed helper — your choice, but it must be COMMITTED and deterministic.)
Marker `r1-regions.done`.
**Testing & validation.** Each chosen window contains ≥1 truth SNV (re-check by
`bcftools view -R win.bed truth.vcf.gz | grep -cv "^#"`); windows within
high-conf BED bounds; README explains counts. No engine runs yet.

---

### TASK R2 — Slice-gate: dual-engine parity on 3×3 regions
**Problem statement.** Prove the whole pipeline (inputs → oracle → ours → compare)
works on real data BEFORE the multi-hour run. Gate for R3.
**Solution theory.** Slice each BAM to the 9 regions (3 HCC1395 + later 3 HG008 may
wait for D3 — HCC1395's 3 suffice to open the gate), run oracle + ours per region,
compare. Byte-parity upstream means any non-MATCH is OUR bug — per §0b stop-and-fix.
**Implementation.** Per region `chr:s-e`:
```bash
bash harness/slice.sh tumor.bam  chr:N-M scratch/g4/slices/tumor.rN.bam
bash harness/slice.sh normal.bam chr:N-M scratch/g4/slices/normal.rN.bam
printf '%s\t%s\t%s\t%s\n' CHR START END GENE > scratch/g4/slices/rN.bed   # match harness BED dialect (see harness headers!)
bash harness/run_oracle.sh scratch/g4/ref/ref.fa "scratch/g4/slices/tumor.rN.bam|scratch/g4/slices/normal.rN.bam" scratch/g4/slices/rN.bed scratch/g4/out/java.rN.vcf TUMOR|NORMAL
bash harness/compare.sh scratch/g4/ref/ref.fa "scratch/g4/slices/tumor.rN.bam|scratch/g4/slices/normal.rN.bam" scratch/g4/slices/rN.bed "TUMOR|NORMAL" 2>&1 | tee scratch/g4/logs/r2.rN.log
```
READ `harness/run_oracle.sh` + `compare.sh` HEADERS first and match their argument
dialects exactly (BED column format, `-z` flag, sample naming). Log every invocation
verbatim into `scratch/g4/logs/cmdline.history`.
Marker `r2-slicegate.done` requires **3× MATCH**.
**Testing & validation.** Three `tee`'d logs each end in MATCH; expected-variant
sanity: oracle VCF contains the truth SNV position(s) for that window
(`grep -F POS java.rN.vcf`). If a window yields zero calls in BOTH engines, swap to
next-ranked window from R1 (document swap in r1/README.md) — empty-vs-empty MATCH is
weak evidence.
**On MISMATCH.** Minimize (smallest failing region), capture diffs, STOP, report per
§0b. Do NOT proceed to R3. Do NOT touch comparators.

---

### TASK R3 — Full-WES dual-engine run (T167 headline)
**Problem statement.** The production-scale credential: whole-exome HCC1395
tumor|normal through both engines, byte-comparable outputs, with wall-time and peak-RSS
recorded.
**Solution theory.** Same commands as R2 with the full-exome high-conf BED as region
input; measure with `/usr/bin/time -v`; engines run sequentially (12 c each wants the
whole box); resumable logging.
**Implementation.**
```bash
REF=scratch/g4/ref/ref.fa
PAIR="scratch/g4/hcc1395/tumor.bam|scratch/g4/hcc1395/normal.bam"
/usr/bin/time -v bash harness/run_oracle.sh $REF "$PAIR" scratch/g4/hcc1395/highconf.bed scratch/g4/out/java.fullwes.vcf "TUMOR|NORMAL" 2> scratch/g4/logs/time.java.txt
/usr/bin/time -v env LC_ALL=C cargo run --release -p vardict-cli -- -G $REF -b "$PAIR" -N "TUMOR|NORMAL" <full-exome flags mirroring run_oracle.sh dialect> > scratch/g4/out/ours.fullwes.var 2> scratch/g4/logs/time.ours.txt
# final compare via compare.sh (it encapsulates the ours-side invocation):
/usr/bin/time -v bash harness/compare.sh $REF "$PAIR" scratch/g4/hcc1395/highconf.bed "TUMOR|NORMAL" 2> scratch/g4/logs/time.compare.txt | tee scratch/g4/logs/r3.result
```
Prefer driving BOTH sides through the harness scripts (single dialect source of
truth); the raw cargo-run line above is documentation of what compare does — use
whichever the harness supports, don't invent flags.
Marker `r3-fullwes.done` on MATCH + captured timings.
**Testing & validation.** `grep MATCH scratch/g4/logs/r3.result`; extract
`Maximum resident set size` + wall clock from both time-files into
`scratch/g4/logs/r3.summary.md`; variant-count sanity: java VCF line count within 2×
of ours (they should be near-identical; gross divergence = investigate before trusting
MATCH — e.g., wrong region dialect silently compared nothing).
**Failure handling.** OOM/kill → reduce concurrent load (confirm deluge/SAB resumed
traffic isn't spiking memory via `free -g`), retry once; repeated failure → report.
Non-MATCH → same discipline as R2.

---

### TASK R4 — HG008 slice validation (T168)
Same pattern as R2 using D3's BAM(s) + HG008 truth-derived regions (from R1-style
selection on the HG008 truth VCF). 3 slices × MATCH. Marker `r4-hg008.done`. Runs
anytime after D3; may interleave with R3 (slice runs are short).

### TASK R5 — Truth-tier concordance (T169)
**Problem statement.** Secondary confirmation that parity-to-Java isn't parity-to-wrong:
score both engines' calls against truth; their precision/recall must be ≈equal (since
outputs match) and reported honestly in absolute terms.
**Solution theory.** `rtg vcfeval --squash-ploidy` (somatic calls vs diploid truth)
for small variants; `truvari bench` for SVs incl. BNDs; restricted to high-conf BED.
Compare OURS-vs-truth vs JAVA-vs-truth: identical inputs ⇒ near-identical PR tables;
any meaningful delta = investigation trigger.
**Implementation.**
```bash
RTG="conda run -n g4-tools rtg"; TRU="conda run -n g4-tools truvari"
$RTG format -o scratch/g4/rtg.sdf scratch/g4/ref/ref.fa
for side in java ours; do
  $RTG vcfeval -b scratch/g4/truth/hcc1395.truth.vcf.gz -c scratch/g4/out/$side.fullwes.vcf(.gz) \
       -t scratch/g4/rtg.sdf --squash-ploidy -e scratch/g4/truth/highconf.bed --all-records \
       -o scratch/g4/eval/$side 2>&1 | tee scratch/g4/logs/vcfeval.$side.log
done
paste <(grep -E "True-pos|False-pos|Precision|Recall|F-score" scratch/g4/eval/java/summary.txt) \
      <(grep -E "True-pos|False-pos|Precision|Recall|F-score" scratch/g4/eval/ours/summary.txt) \
      > scratch/g4/logs/r5.pr-table.txt
$TRU bench -b truth.vcf -c ours.vcf -f ref.fa --includebed highconf.bed -o scratch/g4/bench/ours 2>&1 | tail -20
```
(Adapt filenames: vcfeval wants bgzipped+indexed VCFs — `bgzip -i` + `bcftools index -t`.)
Markers `r5-truth.done`.
**Testing & validation.** PR table exists; per-metric |ours−java| == 0 within rounding
(any nonzero delta → DIVERGENCES entry + STOP for direction); absolute numbers written
into deliverables WITHOUT editorializing ("precision X% vs truth" stated plainly,
Java's number beside ours). Truvari summary saved alongside.

### TASK R6 — CRAM twin piggyback (feeds T156 evidence)
**Problem statement.** §11 claims CRAM support (#249/#380); prove it on real data.
**Solution theory.** CRAM is a lossless-with-reference transform of a validated input;
output must equal the BAM run's output.
**Implementation.**
```bash
conda run -n g4-tools samtools view -C -T scratch/g4/ref/ref.fa scratch/g4/slices/tumor.r2.bam -o scratch/g4/slices/tumor.r2.cram
bash harness/compare.sh scratch/g4/ref/ref.fa "scratch/g4/slices/tumor.r2.cram|scratch/g4/slices/normal.r2.bam" scratch/g4/slices/r2.bed "TUMOR|NORMAL" | tee scratch/g4/logs/r6.log
```
Marker `r6-cram.done` on MATCH (log archived). If compare.sh rejects CRAM inputs,
that itself is finding #T156-partial → report honestly, do not hack around.

---

## DELIVERABLE TASKS

### TASK F1 — Benchmarks & manifest finalization
Copy `r3.summary.md` numbers into `benchmarks/README.md` (new "Full external datasets"
section: dataset, command, wall time, peak RSS, result); ensure `benchmarks/MANIFEST`
covers EVERY downloaded artifact (URL/md5/date/command). Commit both (no data files!).

### TASK F2 — Docs & ledger
- README Validation status: replace "planned, not yet executed" line with actual
  results table (T167/T168 MATCH counts, T169 PR table pointer, R6 CRAM result).
- `docs/upstream-bugs.md`: add entries only for defects surfaced that reproduce with
  the ORACLE alone (manual verification done, non-AI test case attached).
- Ledger append (`.superpowers/sdd/2026-08-21-parity-to-v1.0/progress.md`): per-task
  records + evidence paths.
- Final gates: fmt/clippy/test green; corpus 91/1 unchanged.

---

## DISPATCH GUIDE (for subagent-driven development)
1. **Controller** executes OS-1 personally (privileged host surgery; no sandbox).
2. Dispatch implementers one task at a time in dependency order:
   `D2,T1` (parallel, during drain) → `D1` → `OS-2` → `R1` → `R2` → `R3` ∥ `D3,R4,R6` → `R5` → `F1,F2` → `OS-3`.
   Wait: D3 may run parallel with R3; R4 after D3; R5 strictly after R3+R4.
3. Each dispatch: prompt = "Execute TASK <ID> from
   `docs/plans/2026-08-24-g4-execution.md` VERBATIM; read CLAUDE.md first; check state
   markers; write evidence under scratch/g4/logs/; return commits + gate outputs."
4. Reviewer after each code-touching task (scripts/helpers): spec-adherence +
   secrets-scan + idempotency; reviewers NEVER re-run multi-hour gates.
5. Controller owns milestone transitions OS-2/OS-3 and all commit decisions.
6. Commits: scripts/docs/MANIFEST only. Data, logs, eval outputs stay in scratch/.

## Failure playbook (global)
| Symptom | Action |
|---|---|
| compare.sh MISMATCH anywhere | STOP ladder; minimize; §0b report; fix-forward only via port-loop discipline |
| Download stall/corruption ×2 | Report URL + md5 evidence; try alternate mirror from same index; else abort task |
| Disk <100 GB free | Pause downloads (`mode=pause`), purge nothing, report |
| Container unexpectedly up mid-window | `docker ps` audit + re-stop; check whether a gated cron leaked (grep crontab) |
| QoS daemon found active mid-window | `systemctl --user mask --now` again; investigate what unmasked it before proceeding |
