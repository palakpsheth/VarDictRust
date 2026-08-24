# Tasks: VarDictRust Port (Unified)

> NOTE (2026-08-21): checkboxes below are stale. Authoritative status: DIVERGENCES.md + git log + docs/plans/2026-08-21-parity-to-v1.0.md.

**Input**: Design documents from `/specs/001-vardict-rust-port/`

**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project initialization and basic workspace structure (T001–T014)

- [X] T001 Git branch checkout for init-scaffold
- [X] T002 Create workspace directory layout
- [X] T003 Initialize Cargo workspace manifests in Cargo.toml
- [X] T004 Pin dependencies version specifications in Cargo.toml
- [X] T005 [P] Create core library module stubs in crates/vardict-core/src/lib.rs
- [X] T006 [P] Implement core types and unit tests in crates/vardict-core/src/types.rs
- [X] T007 [P] Create CLI binary stubs in crates/vardict-cli/src/main.rs and cli.rs
- [X] T008 [P] Link CLAUDE.md guidelines file to AGENTS.md
- [X] T009 [P] Create LICENSE and README.md attribution documents
- [X] T010 [P] Create comparative test harness scripts in harness/compare.sh
- [X] T011 Set up Oracle conda environment via harness/setup_oracle.sh
- [X] T012 Seed golden oracle expected outputs in golden/
- [X] T013 [P] Create CI workflow configuration yaml files in .github/workflows/
- [X] T014 [P] Create subagents commands markdown templates in .claude/commands/

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core mathematical, statistical, and data mapping infrastructure (T015–T024)

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

- [X] T015 Implement config structures from CLI arguments in crates/vardict-core/src/config.rs
- [X] T016 [P] Parse BED 3/4-column regions in crates/vardict-core/src/regions.rs
- [X] T017 [P] Parse BED 8-column amplicon regions in crates/vardict-core/src/regions.rs
- [X] T018 Parse region CLI flag option in crates/vardict-core/src/regions.rs
- [X] T019 [P] Implement region tiling window helper in crates/vardict-core/src/regions.rs
- [X] T020 [P] Implement FASTA reader fetching indexed reference bases in crates/vardict-core/src/reference.rs
- [X] T021 [P] Build reference window seed map in crates/vardict-core/src/reference.rs
- [X] T022 Implement indexed BAM reader querying chromosomal ranges in crates/vardict-core/src/bam.rs
- [X] T023 Implement SAM flag, MAPQ, and downsampling filters in crates/vardict-core/src/bam.rs
- [X] T024 Implement CRAM file format decoding utilizing Noodles in crates/vardict-core/src/bam.rs

**Checkpoint**: Foundation ready - user story implementation can now begin in parallel

---

## Phase 3: User Story 1 - Single-Sample (Germline) Variant Calling (Priority: P1) 🎯 MVP

**Goal**: Process a single BAM/CRAM file to call germline variants, outputting VCF matching the oracle's Simple mode (36/38 columns) (T025–T108).

**Independent Test**: Execute `vardictrust` in germline mode against the single-sample test cases and compare outputs using `harness/compare.sh`.

### Tiling & Concurrency (T025–T032)
- [X] T025 [US1] Implement parallel tile runner with ordered output stream in crates/vardict-core/src/runner.rs
- [X] T026 [US1] Implement stable tie-break sorting for deterministic output in crates/vardict-core/src/runner.rs
- [X] T027 [US1] Add CLI flag mapping for filtering arguments (-x, -X, -3, -B) in crates/vardict-core/src/config.rs
- [X] T028 [US1] Add CLI flag mapping for alignment parameters (-o, -O, -I, -z) in crates/vardict-core/src/config.rs
- [X] T029 [US1] Add CLI flag mapping for frequency filtering (-F, -w, -W, -A, -Y, -Z) in crates/vardict-core/src/config.rs
- [X] T030 [US1] Add CLI flag mapping for chimeric/duprate filters (--chimeric, --deldupvar) in crates/vardict-core/src/config.rs
- [X] T031 [US1] Add CLI flag mapping for depth and noise parameters (-mfreq, -nmfreq, -J, -j) in crates/vardict-core/src/config.rs
- [X] T032 [US1] Add CLI flag mapping for standard calling parameters (-i, -p, -DP) in crates/vardict-core/src/config.rs

### CIGAR Modify & Pileup (T033–T058)
- [X] T033 [US1] Implement Cigar type structure and parsing in crates/vardict-core/src/cigar/mod.rs
- [X] T034 [US1] Implement leading/trailing deletion stripping in crates/vardict-core/src/cigar/modify.rs
- [X] T035 [US1] Convert leading/trailing insertions to softclips in crates/vardict-core/src/cigar/modify.rs
- [X] T036 [US1] Implement chimeric softclip removal using seed mapping in crates/vardict-core/src/cigar/modify.rs
- [X] T037 [US1] Implement fixed-point loop for indel cleanups in crates/vardict-core/src/cigar/modify.rs
- [X] T038 [US1] Implement complex fold modifications (combineToCloseToOne/Correct) in crates/vardict-core/src/cigar/modify.rs
- [X] T039 [US1] Implement mis-softclip mismatch recovery in crates/vardict-core/src/cigar/modify.rs
- [X] T040 [US1] Implement leading clip merging (combine_dig_s_dig_m) in crates/vardict-core/src/cigar/modify.rs
- [X] T041 [US1] Verify crash-fixture inputs without panics in crates/vardict-core/src/cigar/modify.rs
- [X] T042 [US1] Implement Variation count trackers in crates/vardict-core/src/variations/mod.rs
- [X] T043 [US1] Define Sclip, Cluster, and Mate tracker maps in crates/vardict-core/src/variations/mod.rs
- [X] T044 [US1] Implement CigarParser pileup entry point scaffold in crates/vardict-core/src/cigar/parse.rs
- [X] T045 [US1] Implement read mismatch gate filter in crates/vardict-core/src/cigar/parse.rs
- [X] T046 [US1] Integrate CigarModifier and softclip length checks in crates/vardict-core/src/cigar/parse.rs
- [X] T047 [US1] Implement match/mismatch base loops in crates/vardict-core/src/cigar/parse.rs
- [X] T048 [US1] Implement multi-nucleotide variant (MNV) growth check in crates/vardict-core/src/cigar/parse.rs
- [X] T049 [US1] Implement deletion (D) event counting in crates/vardict-core/src/cigar/parse.rs
- [X] T050 [US1] Implement insertion (I) event counting and homopolymer shift check in crates/vardict-core/src/cigar/parse.rs
- [X] T051 [US1] Implement splice (N) event tracking in crates/vardict-core/src/cigar/parse.rs
- [X] T052 [US1] Implement softclip read mapping extraction in crates/vardict-core/src/cigar/parse.rs
- [ ] T053 [US1] Test simple mode pre-realignment intermediate output for chromosome 10 in crates/vardict-core/src/cigar/parse.rs
- [ ] T054 [US1] Verify depth and strand bias column parity in raw simple tables in crates/vardict-core/src/cigar/parse.rs
- [ ] T055 [US1] Test insertion/deletion variant counts parity before realignment in crates/vardict-core/src/cigar/parse.rs
- [ ] T056 [US1] Verify homopolymer alt base filtering parity in crates/vardict-core/src/cigar/parse.rs
- [ ] T057 [US1] Test mate and soft-clipped mapping metrics parity in crates/vardict-core/src/cigar/parse.rs
- [ ] T058 [US1] Run Tier-1 pre-realignment differential comparisons using harness/compare_var.py

### Realignment (T059–T075)
- [X] T059 [US1] Implement consensus check and correct count utilities in crates/vardict-core/src/variations/util.rs
- [X] T060 [US1] Implement ismatch base mismatch filter checks in crates/vardict-core/src/realign.rs
- [X] T061 [US1] Implement breakpoint scanning loops (findbp) in crates/vardict-core/src/realign.rs
- [X] T062 [US1] Implement flanking mismatch checks (findMM3) in crates/vardict-core/src/realign.rs
- [X] T063 [US1] Implement deletion realignment logic in crates/vardict-core/src/realign.rs
- [X] T064 [US1] Implement insertion realignment logic in crates/vardict-core/src/realign.rs
- [X] T065 [US1] Implement large-deletion softclip realignment using sub-pipelines in crates/vardict-core/src/realign.rs
- [X] T066 [US1] Implement large-insertion softclip realignment and duplication tracking in crates/vardict-core/src/realign.rs
- [X] T067 [US1] Define shared structural variant markers in crates/vardict-core/src/sv.rs
- [X] T068 [US1] Integrate deletion realignment loop and test EGFR deletions in crates/vardict-core/src/realign.rs
- [X] T069 [US1] Integrate insertion realignment loop and test EGFR insertions in crates/vardict-core/src/realign.rs
- [ ] T070 [US1] Test softclip realignment for candidate deletion records in crates/vardict-core/src/realign.rs
- [ ] T071 [US1] Test softclip realignment for candidate insertion records in crates/vardict-core/src/realign.rs
- [ ] T072 [US1] Validate realigned insertion duplication count parity in crates/vardict-core/src/realign.rs
- [ ] T073 [US1] Verify structural variant breakpoint coordinates in crates/vardict-core/src/realign.rs
- [ ] T074 [US1] Run Tier-1 full realignment comparisons using harness/compare_var.py
- [ ] T075 [US1] Ensure all realignment heuristics match the oracle on testdata cases

### Finalization (T076–T088)
- [X] T076 [US1] Define Variant output structure fields in crates/vardict-core/src/variations/mod.rs (merged into mod.rs)
- [X] T077 [US1] Implement createVariant and alternate frequency checks in crates/vardict-core/src/tovars.rs
- [X] T078 [US1] Implement stable sorting of variants using quality/coverage keys in crates/vardict-core/src/tovars.rs
- [X] T079 [US1] Implement MSI motif scanner and shift3 coordinates in crates/vardict-core/src/tovars.rs
- [X] T080 [US1] Parse variant allele description grammar in crates/vardict-core/src/tovars.rs
- [X] T081 [US1] Implement isGoodVar threshold and quality metrics filters in crates/vardict-core/src/variations/variant.rs
- [X] T082 [US1] Implement SimplePostProcess caller module in crates/vardict-core/src/postprocess/simple.rs
- [X] T083 [US1] Implement 36-column printer logic for Simple mode in crates/vardict-core/src/output/table.rs
- [X] T084 [US1] Test 36-column output columns layout and order in crates/vardict-core/src/output/table.rs
- [X] T085 [US1] Verify allele count and frequency columns parity in crates/vardict-core/src/output/table.rs
- [X] T086 [US1] Test MSI and shift3 columns output parity in crates/vardict-core/src/output/table.rs
- [X] T087 [US1] Validate mean quality and mapping quality columns parity in crates/vardict-core/src/output/table.rs
- [X] T088 [US1] Run Tier-1 germline differential comparisons using harness/compare_var.py

### Statistical Parity (Fisher) (T089–T096)
- [X] T089 [US1] Implement hypergeometric probability density function in crates/vardict-core/src/stats/fisher.rs
- [X] T090 [US1] Implement log-sum-exp hypergeometric cumulative calculations in crates/vardict-core/src/stats/fisher.rs
- [X] T091 [US1] Implement two-sided Fisher's Exact test p-value calculations in crates/vardict-core/src/stats/fisher.rs
- [X] T092 [US1] Implement one-sided Fisher p-value helper bounds in crates/vardict-core/src/stats/fisher.rs
- [X] T093 [US1] Implement Brent's root-finding solver (zeroin) in crates/vardict-core/src/stats/zeroin.rs
- [X] T094 [US1] Implement conditional Maximum Likelihood Estimation of odds ratios in crates/vardict-core/src/stats/fisher.rs
- [X] T095 [US1] Implement statistical cross-checks against fishers_exact crate in crates/vardict-core/src/stats/fisher.rs
- [X] T096 [US1] Implement 38-column output layout and append Fisher calculations in crates/vardict-core/src/output/table.rs

### VCF Output Validation (T097–T108)
- [X] T097 [US1] Format VCF v4.2 headers with contig lengths in crates/vardict-core/src/output/vcf.rs
- [X] T098 [US1] Port quality filtering logic (d3, v2, f0.02) in crates/vardict-core/src/output/vcf.rs
- [X] T099 [US1] Port mapping quality filtering logic (q22.5, Q10) in crates/vardict-core/src/output/vcf.rs
- [X] T100 [US1] Port strand bias filter checks (Bias) in crates/vardict-core/src/output/vcf.rs
- [X] T101 [US1] Port mismatch and gap filter rules (NM5.25, p8) in crates/vardict-core/src/output/vcf.rs
- [X] T102 [US1] Port position standard deviation rules (pSTD) in crates/vardict-core/src/output/vcf.rs
- [X] T103 [US1] Port homopolymer and MSI filter checks (MSI12) in crates/vardict-core/src/output/vcf.rs
- [X] T104 [US1] Port signal-to-noise filter thresholds (SN1.5) in crates/vardict-core/src/output/vcf.rs
- [X] T105 [US1] Implement VCF INFO and FORMAT record serialization in crates/vardict-core/src/output/vcf.rs
- [X] T106 [US1] Wire CLI coordinator to simple variant calling pipeline in crates/vardict-cli/src/main.rs
- [X] T107 [US1] Run Tier-2 simple mode differential testing using harness/compare.sh
- [X] T108 [US1] Snapshot expected simple VCF outputs in golden/

**Checkpoint**: At this point, User Story 1 is fully functional and testable independently.

---

## Phase 4: User Story 2 - Somatic (Paired Tumor/Normal) Variant Calling (Priority: P2)

**Goal**: Process paired alignments and output somatic status classifications matching the oracle's Somatic mode (55/61 columns) (T109–T124).

**Independent Test**: Execute `vardictrust` in somatic mode (`-b "tumor.bam|normal.bam"`) and compare outputs using `harness/compare.sh`.

- [ ] T109 [US2] Implement two-track somatic alignment stream loader in crates/vardict-core/src/modes/somatic.rs
- [ ] T110 [US2] Implement 55-column printer output layout for Somatic mode in crates/vardict-core/src/output/table.rs
- [ ] T111 [US2] Implement 61-column Fisher-included Somatic output table in crates/vardict-core/src/output/table.rs
- [ ] T112 [US2] Implement read noise filtering (isNoise) in crates/vardict-core/src/variations/mod.rs
- [ ] T113 [US2] Implement determinateType somatic classification rules in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T114 [US2] Implement combineAnalysis joint re-pileup workflow in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T115 [US2] Implement sample acceptance dispatch logic in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T116 [US2] Implement somatic 3-Fisher calculations matching R tests in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T117 [US2] Port var2vcf_paired filtering and INFO/FORMAT tags to VCF writer in crates/vardict-core/src/output/vcf.rs
- [ ] T118 [US2] Validate tumor-only calling somatic output columns parity in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T119 [US2] Validate normal-only calling somatic output columns parity in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T120 [US2] Test joint tumor-normal somatic classification rules parity in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T121 [US2] Test LikelySomatic and LikelyLOH filter thresholds parity in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T122 [US2] Test somatic p-value and odds ratio float values parity in crates/vardict-core/src/postprocess/somatic.rs
- [ ] T123 [US2] Run Tier-2 somatic paired mode differential testing using harness/compare.sh
- [ ] T124 [US2] Snapshot expected somatic paired VCF outputs in golden/

**Checkpoint**: At this point, User Stories 1 AND 2 should both work independently.

---

## Phase 5: User Story 3 - Amplicon Mode and Splicing/SV Support (Priority: P3)

**Goal**: Process PCR-based amplicon sequencing datasets or detect structural and splicing variants matching oracle output format (T125–T146).

**Independent Test**: Execute `vardictrust` in amplicon/splicing mode and compare outputs using `harness/compare.sh`.

### Amplicon Mode (T125–T134)
- [ ] T125 [US3] Implement amplicon overlap alignment soft-filtering in crates/vardict-core/src/cigar/parse.rs
- [ ] T126 [US3] Implement AmpliconMode coordinate-to-amplicon mapping in crates/vardict-core/src/pipeline/amplicon.rs
- [ ] T127 [US3] Implement AmpliconPostProcess bias flag evaluation in crates/vardict-core/src/postprocess/amplicon.rs
- [ ] T128 [US3] Implement 38/40-column amplicon printer output layouts in crates/vardict-core/src/output/table.rs
- [ ] T129 [US3] Test amplicon overlapping read clipping heuristics parity in crates/vardict-core/src/postprocess/amplicon.rs
- [ ] T130 [US3] Test amplicon bias indicators and target coverage metrics parity in crates/vardict-core/src/postprocess/amplicon.rs
- [ ] T131 [US3] Test 38-column and 40-column amplicon intermediate outputs parity in crates/vardict-core/src/output/table.rs
- [ ] T132 [US3] Validate amplicon mode variant calling against oracle using harness/compare.sh
- [ ] T133 [US3] Snapshot expected amplicon VCF outputs in golden/
- [ ] T134 [US3] Ensure amplicon primer coordinates match BED regions during parsing in crates/vardict-core/src/pipeline/amplicon.rs

### Structural Variants & Splicing (T135–T146)
- [x] T135 [US3] Define Structural Variant tracking data structures in crates/vardict-core/src/sv/mod.rs
- [x] T136 [US3] Implement DEL/INV structural variant search processing in crates/vardict-core/src/sv/mod.rs
- [ ] T137 [US3] Implement SplicingMode splice count tracking in crates/vardict-core/src/pipeline/splicing.rs
- [ ] T138 [US3] Add SV-specific INFO tags (SVTYPE, SVLEN) to VCF serializer in crates/vardict-core/src/output/vcf.rs
- [x] T139 [US3] Validate SV candidate deletion structural mappings parity in crates/vardict-core/src/sv/mod.rs
- [ ] T140 [US3] Validate SV candidate inversion structural mappings parity in crates/vardict-core/src/sv/mod.rs
- [ ] T141 [US3] Validate split-read and span-pair alignment counts parity in crates/vardict-core/src/sv/mod.rs
- [ ] T142 [US3] Test SV duplication rate calculations parity in crates/vardict-core/src/sv/mod.rs
- [ ] T143 [US3] Test splicing variant splice counts output parity in crates/vardict-core/src/pipeline/splicing.rs
- [ ] T144 [US3] Run Tier-2 structural variant differential testing using harness/compare.sh
- [ ] T145 [US3] Run truvari benchmarking comparison for SV structural metrics
- [ ] T146 [US3] Snapshot expected SV/splicing VCF outputs in golden/

**Checkpoint**: All user stories should now be independently functional.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Performance optimization, robustness checks, and release management (T147–T174)

### Performance & Scaling (T147–T154)
- [ ] T147 Set up performance benchmarking suite utilizing criterion in benches/engine.rs
- [ ] T148 Feature-gate AVX2-accelerated block-aligner implementation in crates/vardict-core/src/realign.rs
- [ ] T149 Optimize memory allocations reusing tile buffers in crates/vardict-core/src/runner.rs
- [ ] T150 Benchmark peak RSS and tune tiling coordinates sizing boundaries
- [ ] T151 Test parallel thread throughput scaling metrics
- [ ] T152 Verify thread concurrency safety on multi-thread tile runners
- [ ] T153 Profile CPU hotspots under multi-threaded execution
- [ ] T154 Document execution scaling benchmarks in benchmarks/README.md

### Robustness & Systems (T155–T166)
- [ ] T155 Run tests against the robustness crash corpus in tests/robustness.rs
- [ ] T156 Verify end-to-end CRAM index and alignment reading parity using harness/compare.sh
- [ ] T157 Enable CI workflows checks for cargo formatting and audit
- [ ] T158 Configure cargo clippy warnings rules in github actions yaml
- [ ] T159 Validate locale-independent float parsing and formatting
- [ ] T160 Validate correct contig lengths assignment from reference fai files
- [ ] T161 Set up cargo deny rules configurations in deny.toml
- [ ] T162 Configure cargo dist release binary compiler configurations in Cargo.toml
- [ ] T163 Set up crates.io trusted publisher GitHub actions workflows
- [ ] T164 Implement release-plz workflows for automated changelog updates
- [ ] T165 Complete contribution guide documents in CONTRIBUTING.md
- [ ] T166 Document supported flag matrix and known issues in README.md

### System-wide Validation & Release (T167–T174)
- [ ] T167 Run manual verification sweep on full SEQC2 WES dataset
- [ ] T168 Run manual verification sweep on pancreatic cancer somatic datasets
- [ ] T169 Validate precision-recall rates comparison against GIAB validation sets
- [ ] T170 Execute full §15 acceptance checks list sweep
- [ ] T171 Verify final README provenance details completeness
- [ ] T172 Track remaining upstream bugs and document test cases
- [ ] T173 Release v1.0.0 cargo packages to crates.io and GitHub
- [ ] T174 Document post-1.0 features roadmap in ROADMAP.md

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately.
- **Foundational (Phase 2)**: Depends on Setup completion - BLOCKS all user stories.
- **User Stories (Phase 3+)**: All depend on Foundational phase completion.
  - Can proceed sequentially in priority order (P1 $\rightarrow$ P2 $\rightarrow$ P3).
- **Polish (Phase 6)**: Depends on all user stories being complete.

---

## Parallel Example: User Story 1

```bash
# Launch parallel unit tests for User Story 1:
Task: "T034 - Implement leading/trailing deletion stripping in crates/vardict-core/src/cigar/modify.rs"
Task: "T035 - Convert leading/trailing insertions to softclips in crates/vardict-core/src/cigar/modify.rs"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001–T014)
2. Complete Phase 2: Foundational (T015–T024)
3. Complete Phase 3: User Story 1 (Germline MVP, T025–T108)
4. **STOP and VALIDATE**: Validate germline output parity against oracle using `harness/compare.sh`.

### Incremental Delivery

1. Setup + Foundation $\rightarrow$ workspace ready.
2. Add User Story 1 $\rightarrow$ Test germline $\rightarrow$ Deploy/Demo (MVP ready).
3. Add User Story 2 $\rightarrow$ Test somatic $\rightarrow$ Deploy/Demo.
4. Add User Story 3 $\rightarrow$ Test amplicon/SVs $\rightarrow$ Deploy/Demo (Full parity reached).
