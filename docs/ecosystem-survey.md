# Ecosystem Survey — AI-assisted Rust rewrites of bioinformatics tools

**Author / Maintainer:** [Palak Sheth](https://www.linkedin.com/in/palak-p-sheth/) <palakpsheth@gmail.com>  
**Repository:** [https://github.com/palakpsheth/VarDictRust](https://github.com/palakpsheth/VarDictRust)  
**Governing Manifesto:** [rewrites.bio](https://rewrites.bio)  

> **Point-in-time snapshot: 2026-08-24.** Compiled from web/GitHub sweeps (sources
> listed at the end; claims reflect pages as fetched that day). This document is
> provenance for decisions consolidated into project specifications and external validation reports in
> [`docs/validation/omnibus-empirical-validation.md`](validation/omnibus-empirical-validation.md).
> Weaknesses are given equal standing to strengths — see §5.

## 1. The wave

AI coding assistants crossed a threshold where a domain expert can reimplement a
mature tool in days. **rewrites.bio** (Seqera, Apr 2026) is the movement's anchor
manifesto. Its principles:

1. **Credit the original authors** — README, papers, downstream reports; citations flow to the original.
2. **Emulate exactly** — byte-for-byte for deterministic tools; scientist-defined float tolerance otherwise; "everything counts": headers, column order, file naming.
3. **Be transparent about AI** — which tools, what role, how correctness was verified, where validation is thin.
4. **Think big / work small** — rethink architecture (fold pipeline steps, kill intermediate I/O), but implement in tight loops with output comparison at every step ("AI agents are unreliable for correctness — validate every step").
5. **Test and benchmark with real data**, documented hardware/datasets/commands.
6. **Build only what you need** — unsupported features must fail loudly pointing at the original; never silently ignore flags.
7. **Pin versions and document** — equivalence claims are about an exact version; document commands, data, dates, upstream-update policy.
8. **Stewardship** — maintain what you ship; preserve hot-swappable compatibility; open source; contribute upstream responsibly (verify manually; *never* AI-generated bug-report evidence).

Origin story: RustQC began as dupRadar rewrite over lunch, became a 15-tool QC suite at >60× speedup; weeks of parity iteration followed the weekend beta ("Stringency in validation will be key to community adoption and trust").

## 2. Exemplar rewrites and their evidence models

| Project | Replaces | Evidence model |
|---|---|---|
| [bwa-mem4](https://github.com/IPNP-BIPN/bwa-mem4) | bwa-mem2 2.3 | Byte-identity acceptance gate from first commit; whole-genome GIAB HG002 byte-identical SAM (353M SE + 707M PE records); 2000-job property test on SW kernel; release binaries must rebuild committed fixtures byte-identically per platform |
| [Bismark Rust suite](https://felixkrueger.github.io/Bismark/rust/overview/) | Perl Bismark (12 tools) | Byte-identical contract v0.25.1; worker-invariant `--multicore`; deliberate documented divergences; opt-in non-identical modes concordance-gated |
| [rammap](https://github.com/jwanglab/rammap) | minimap2 | Identical output; found 2 long-standing minimap2 bugs → fixed upstream v2.31 |
| [rustar-aligner](https://scverse.org/rustar-aligner/) | STAR | Public compatibility scoreboard; 99.88% mate faithfulness labeled "tie-adjusted"; 396 tests |
| [tama-rs](https://github.com/sguizard/tama-rs) | TAMA (Python 2) | Byte-identical goldens; dict-iteration-order exceptions compared as sets and documented |
| [RustQC](https://github.com/seqeralabs/RustQC) | 15 RNA-seq QC tools | nf-test snapshot harness as agent feedback loop; public AI & Provenance page incl. known gaps |
| [fastqc-compliant-rs](https://github.com/henriksson-lab/fastqc-compliant-rs) | Java FastQC | Golden char-for-char files; published parity matrix; refuses replacement claims until matrix green |
| [FastQC-Rust](https://github.com/ewels/FastQC-Rust) | Java FastQC | `UPSTREAM.toml` version pin; nightly CI bot files issues on upstream releases; equivalence report CI |
| [fgumi](https://github.com/fulcrumgenomics/fgumi) | fgbio (Scala) | 25× speedup, 100% output-reproduction claim |
| [fastars](https://github.com/necoli1822/fastars) | fastp | 100% identical sequences verified byte-for-byte; honest "limited samples" warning |

Also notable: sshash-rs/piscem (Patro, C++→Rust in <4 days, 100% concordance);
henriksson-lab's port program (minimap2, BLAST, HMMER, MUSCLE, DIAMOND, Kraken2…);
umi-collapse-rs, dist2mst-rs (~900×), EVidenceModeler_rust, bowtie2-rs ("parity
gates before core aligner work").

## 3. Community feedback

**Pro** (Heng Li, ["The AI Rewrite Dilemma"](https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma)): order-of-magnitude speedups for Python/R tools; static Rust binaries end dependency hell; welcomes rewrites of his own tools. Reversibility makes adoption risk-free — swap in, confirm on one run, roll back if needed.

**Con/skepticism** (Kai Blin thread): *"Who will maintain vibe-rewritten tools? Original maintainers won't."* Unmergeable AI PRs burden maintainers; fragmentation risk — subtle rounding/multimap differences poison cross-study comparability.

**Deeper critique** ([Clay McLeod](https://claymcleod.dev/blog/2026-07-13-beyond-bioinformatics-rewrites)): emulation preserves historical constraints without separating durable science from accidents; rewrites are seductively easy AI projects; noodles (spec-first foundation) may compound more than any single rewrite; **"a tool author can delegate implementation, but cannot delegate understanding"**; wants provenance-as-runtime-artifact.

**Counter-example teardown** ([polars-bio on RastQC](https://biodatageeks.org/polars-bio/blog/2026/07/04/cheap-to-rewrite-expensive-to-get-right-the-fastqc-compatible-trap/)): claimed FastQC-compatible with wrong core outputs, non-thread-invariant sampling, "speedup" that was parallelism-not-port (slower than the JVM single-threaded, silently gated above 50 MB). Thesis: *a rewrite is a validation contract; cheap generation needs disciplined validation* — spec-driven + test-driven development, golden assertions not header greps, partition-invariance checks, like-for-like benchmarks.

**Trust levers that convinced skeptics:** byte-diff against a pinned oracle on real data > any claim; published parity matrices; thread/partition-invariance proofs; manually-verified upstream bug reports; gating artifacts (not just source) before release.

## 4. Scientific Variant Caller Ecosystem Comparison

| Dimension | VarDictRust v1.0.x | VarDictJava v1.8.3 | GATK Mutect2 (v4.5+) | Illumina Strelka2 (v2.9+) | Octopus (v0.7+) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Language / Runtime** | Single static Rust binary (`musl`/`glibc`) | Java 8 + R 3.x + Perl 5 pipeline | Java JVM (GATK framework) | C++11 + Python wrapper | C++14 native |
| **Memory Footprint** | Bounded ($O(1)$ halo tiling, `1.27 GiB` on 249 Mbp `chr1`) | Unbounded ($O(N)$ depth, OOM prone) | High (JVM heap 8–16 GB) | Moderate (4–8 GB) | High (Haplotype graph 8–16 GB) |
| **Execution Model** | Rayon work-stealing, thread-deterministic | Multi-process pipe, non-deterministic | Multi-threaded JVM | Multi-threaded task graphs | Multi-threaded DAG |
| **Calling Modalities** | Germline, Somatic, Amplicon, Splicing/SV | Germline, Somatic, Amplicon, Splicing/SV | Somatic, Mitochondrial | Germline, Somatic | Germline, Somatic, De Novo |
| **Complex Realignment**| In-engine Smith-Waterman + soft-clip assembly | In-engine Smith-Waterman + soft-clip assembly | Local Haplotype assembly | Prefix tree / HMM alignment | Haplotype-based HMM |
| **Low-AF Sensitivity**  | Ultra-high ($AF \ge 0.1\%$, down to 2 reads) | Ultra-high ($AF \ge 0.1\%$, down to 2 reads) | Moderate ($AF \ge 1\%$) | Moderate ($AF \ge 1\%$) | High ($AF \ge 0.5\%$) |
| **Targeted Amplicon**  | Primer trimming (`-a`), amplicon bias filter | Primer trimming (`-a`), amplicon bias filter | Requires external trimming | External trimming | External trimming |
| **Upstream Parity**    | Byte-for-byte on valid inputs ([rewrites.bio](https://rewrites.bio)) | Baseline reference | N/A | N/A | N/A |

### How VarDictRust Meets rewrites.bio Principles

| Field Practice | VarDictRust Implementation |
| :--- | :--- |
| **Emulate exactly, bugs included** | [AGENTS.md §5 Bug-Replication Policy](../AGENTS.md#5-emulate-exactly-the-rewritesbio-prime-directive) and [`DIVERGENCES.md`](../DIVERGENCES.md) ledger. Upstream bugs replicated in v1.0.0, slated for opt-in flags in v1.1.0 per [`ROADMAP.md`](../ROADMAP.md). |
| **Pinned oracle + differential testing** | VarDictJava v1.8.3 @`0886261` conda environment; 92-case corpus; live-oracle VCF sweeps. |
| **Work small, per-function evidence** | [AGENTS.md §7 Port Loop](../AGENTS.md#7-the-port-loop-rewritesbio-work-small) with exact Java source code citations. |
| **Fail loudly on unsupported features** | [AGENTS.md §12 Unsupported Error Policy](../AGENTS.md#12-unsupported-features-error-loudly) implemented in [`cli.rs`](../crates/vardictrust/src/cli.rs) via `UNSUPPORTED_FLAGS`. |
| **Real-data, like-for-like benchmarks** | Documented in [`benchmarks/README.md`](../benchmarks/README.md); strictly avoided parallelism-only speedup claims. |
| **Thread & partition invariance** | Formally verified in [`engine_behavior.rs`](../crates/vardictrust-core/tests/engine_behavior.rs). |
| **Upstream bugs, manually verified** | Documented in [`docs/upstream-bugs.md`](upstream-bugs.md) ([UB-1](upstream-bugs.md#ub-1--p-spurious-all-zero-record-at-position-00) through [UB-7](upstream-bugs.md#ub-7-multi-allelic-site-single-winner-collapsing-var2vcf_validpl135-205)). |
| **Credit & AI transparency** | Comprehensive attribution and disclosures in [`README.md`](../README.md) and [`benchmarks/evidence/`](../benchmarks/evidence/). |

---

## 5. Architectural Hardening & Verification Matrix

| Area | Historical Gap | Engine Hardening Measure | Verification Evidence |
| :--- | :--- | :--- | :--- |
| **Hostile CIGAR Fuzzing** | Crash risks in CIGAR modifier | Bounds-checked CIGAR manipulation in `crates/vardictrust-core/src/cigar/` | `crates/vardictrust/tests/robustness.rs` covering PR #402 corrupt topology |
| **CLI Flag Parity** | Incomplete CLI flag documentation | Explicit 89-flag Diátaxis reference categorized into 7 functional groups (plus [`docs/migration.md`](migration.md)) | Section 4 in [`README.md`](../README.md) |
| **Release Artifact Packaging** | Ungated release binaries | PIE static binary validation gate comparing musl to glibc outputs | [`harness/bin_equiv.sh`](../harness/bin_equiv.sh), [`musl-ci.json`](../benchmarks/evidence/musl-ci.json) |
| **Biological Truth Concordance** | Lack of physical truth-set evaluation | Full multi-cohort RTG `vcfeval` truth concordance suite | 5 cohorts in [`docs/validation/omnibus-empirical-validation.md`](validation/omnibus-empirical-validation.md) |
| **Upstream Drift Defense** | Unpinned upstream dependency risk | Pinned reference oracle to VarDictJava v1.8.3 (`commit 0886261`) | Sealed environment specification in [`harness/setup_oracle.sh`](../harness/setup_oracle.sh) |
| **Bioconda Packaging** | No automated Conda distribution | Canonical in-tree Bioconda recipe under `packaging/bioconda/recipes/vardictrust/` | [`harness/tests/test_bioconda_recipe.py`](../harness/tests/test_bioconda_recipe.py) |

## 6. Sources (as fetched 2026-08-24)

- https://rewrites.bio · https://github.com/seqeralabs/rewrites.bio · https://seqera.io/blog/rewrites-bio/ · https://seqera.io/blog/rustqc/
- https://seqeralabs.github.io/RustQC/about/ai-statement/
- https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma · https://lh3.github.io/2026/07/04/minibwa-is-the-new-bwa
- https://claymcleod.dev/blog/2026-07-13-beyond-bioinformatics-rewrites
- https://biodatageeks.org/polars-bio/blog/2026/07/04/cheap-to-rewrite-expensive-to-get-right-the-fastqc-compatible-trap/
- https://scverse.org/rustar-aligner/ · https://github.com/IPNP-BIPN/bwa-mem4/releases/tag/v4.0.0
- https://felixkrueger.github.io/Bismark/rust/overview/ · https://github.com/henriksson-lab/fastqc-compliant-rs
- https://github.com/ewels/FastQC-Rust · https://github.com/sguizard/tama-rs · https://github.com/necoli1822/fastars · https://github.com/jwanglab/rammap
- Phil Ewels, Nextflow Summit 2026 talk + launch threads (Kai Blin rebuttals on Bluesky/mu.social)

## 7. Addendum 2026-10-11: vardictcpp third arm (this section only is newer than the snapshot)

Every Java-vs-Rust comparison in this repo now runs beside **vardictcpp 3
`hc52dbad_0` bioconda** (same `vardict-oracle` conda env as the VarDictJava
oracle) as an **advisory third arm**: `CPP: MATCH/MISMATCH/SKIP` lines and
`compare_cpp.txt` artifacts, evidence `extra` keys (`cpp_version`, `cpp_head`,
`cpp_wall_s`, `cpp_maxrss_kb`, `cpp_f1`, `cpp_skip_reason`), `CPP_SKIP=1`
escape. The release gate stays java==rust; cpp divergences are normalized in
[`docs/upstream-bugs.md`](upstream-bugs.md#cpp-divergences),
not silently dropped. Galaxy wrapper distribution of VarDictRust itself stays a
follow-up (see [`CHANGELOG.md`](../CHANGELOG.md)).

---

<a id="critical-review-v107"></a>
## 8. Addendum 2026-11-08: Clinical, Scientific & Engineering Critical Review (v1.0.7)

This section synthesizes a 5-domain critical engineering and clinical oncology survey across:
1. **Fulcrum Genomics & `rewrites.bio` Rust Rewrites** (`fgumi`, `perbase`, `ferro-hgvs`, `riker`, `bwa-mem4`, `Bismark`, `FastQC-Rust`, `RustQC`, `vardictcpp`)
2. **PacBio Native Rust Suite & Greenfield Rust Callers** (`sawfish`, `trgt`, `HiPhase`, `paraphase`, `pb-CpG-tools`, `aardvark`, `fibertools-rs`, `varlociraptor`, `kanpig`, `modkit`, `oarfish`, `piscem`)
3. **Peer-Reviewed Clinical Oncology & Variant Calling Literature (2016–2026)**
4. **Community Discourse Across LinkedIn, Reddit (`r/bioinformatics`, `r/rust`), Biostars, and GitHub Issues/PRs**
5. **Internal Line-by-Line Codebase & Performance Audit of `VarDictRust` v1.0.7**

### 8.1 SWOT Analysis (`VarDictRust` v1.0.7)

| Quadrant | Core Findings |
| :--- | :--- |
| **Strengths** | • **100% Oracle Parity:** Bitwise/float ($\varepsilon=10^{-6}$) equivalence to `VarDictJava v1.8.3` across 92 Tier-1 fixtures, 92/92 full-chromosome validation runs (70 Tier-2 full-VCF comparisons + 22 Tier-3 count/PASS-site/`vcfeval` checks; [`chr-parity-v107.md`](validation/chr-parity-v107.md)), 17 modalities, 18 multi-platform slices, and 5 NIST GIAB / FDA SEQC2 truth sets.<br/>• **Bounded $O(1)$ Memory via Halo Tiling:** Flat **1.11× RSS scaling** on `chr1` (`1.15 GiB` at 10 Mbp $\to$ `1.27 GiB` at 249 Mbp vs Java OOM and `vardictcpp` `24.34 GiB` unchunked RSS); only VarDict engine completing 364× HG008 paired WGS (`329 GB` BAM) in `1h45m` (`4.99–12.43 GiB` RSS).<br/>• **Memory-Safe Core & HTSlib-Free I/O Stack:** `#![forbid(unsafe_code)]` in `vardictrust-core` + HTSlib-free `noodles 0.118` I/O with native `noodles-bgzf 0.52` (`libdeflate` C backend + `zlib-rs`), eliminating all `modifyCigar` index crashes and `hts-sys` bindgen breakage.<br/>• **Single-Binary Consolidation & Supply Chain:** Fuses `Java -> R -> Perl` into one static binary (`cargo-dist` SLSA L3, Crates.io, Bioconda, Seqera Wave, `nf-core/modules` `#13343`).<br/>• **Clinical Niche:** Best-in-class sensitivity for complex/composite indels and MNVs (*BRCA1/2*, *EGFR* Ex19del/Ex20ins, *FLT3*-ITD), native 8-column BED amplicon primer/overlap bias (`-a` / `AMPBIAS`), and unified Somatic + Germline + LOH classification. |
| **Weaknesses** | • **Resolved `v1.0.7` Hotspots (`PR #98`–`#107` on `main` ahead of `v1.0.8`):** Multi-tile `scan_multi_tile_max_read_lengths` pre-scans are now parallelized across `rayon`; `modify_cigar` now fast-paths single-`M` and pure `S`/`M` reads in-place without `String`/regex round-trips; `Reference::get_base` inlines primary-region lookup; `BGZF_BLOCK_CACHE` expanded from 8 to 32 shards; `FisherExact` eliminated `support: Vec<i32>` and `sum_elements` heap allocations; `cli.rs` `--version` license text matches `MIT AND GPL-2.0-or-later`; and native `--output .vcf.gz` + `--tbi`/`--csi` + all 17 `var2vcf` `--vcf-*` flags + `--completions`/`--manpage` landed.<br/>• **Remaining Engine & Scaling Gaps:** `DEFAULT_UNTILED_SPAN_MAX = 1_000_000` still leaves 100 kbp–1 Mbp BED intervals untiled by default; `BamFeed` materializes `BamRecord` before coordinate/flag pre-filtering; `sv.rs` `find_inv_sub` still runs eager mate-window `partial_pipeline` (`sv-lazy-mate-pileup`, deferred to `v1.1.0` `--modern` `#95`).<br/>• **GPL-2.0-or-later Contamination:** 5 small R `nmath` files (~350 LoC in `stats/rmath/`) force `MIT AND GPL-2.0-or-later` workspace licensing until cleanroom replacement in `v1.1.0` (`#92`).<br/>• **Inherited `v1.8.3` UX & Clinical Quirks:** Default 13-column UCSC `refGene` settings (`-c 3 -S 7 -E 8 -s 10 -e 11 -g 13`, with 4-column standard BED and 8-column amplicon BED auto-detected) require explicit `-c 1 -S 2 -E 3 -z` on 3-, 5-, and 6-column UCSC BEDs (now emitting an actionable `VardictError::Parse` hint); `-A` flag collision (`sv_min_af` vs `var2vcf_valid.pl -A` all-alleles `UB-7`); `VD > DP` (`#92`); elevated raw low-VAF FDR without UMI/duplex or `F1R2`/`F2R1` orientation bias filtering. |
| **Opportunities** | • **Surpass `vardictcpp` on All Germline Workloads (`v1.0.x` Tier 1):** Output-neutral fast-paths (`modify_cigar` single-`M` fast-path preserving endpoint mismatch clipping, parallel `max_read_length` scan, 100 kbp tiling threshold, cached `&LoadedRegion`, 32-shard BGZF cache) will cut germline wall times by 30–45%.<br/>• **Restore Pure `MIT` License (`v1.1.0`):** Replace `stats/rmath/` (~350 LoC) with a cleanroom Loader (2000) / Apache-2.0 saddle-point implementation.<br/>• **Fulcrum & PacBio Engineering Convergence:** Native `-o out.vcf.gz` + `.tbi` via `noodles::bgzf::io::MultithreadedWriter` (`aardvark` pattern), undecoded `noodles::bam::Record` pre-filtering (`fgumi-raw-bam` / `oarfish` pattern), and `--stats-json` telemetry (`sawfish` / `fgumi` pattern).<br/>• **Opt-In `--fix` & `--preset` Evolution (`v1.1+` / `v2.0`):** Ship `--fix vd-dp-reconcile`, `--fix ub7` / `--all-alleles`, `--fix 378` (paired somatic amplicon), `--fix eqx-read-length`, `--fix sv-lazy-mate-pileup`, plus streaming UMI/duplex consensus (`--umi-tag RX --duplex`), `F1R2`/`F2R1` FFPE filter, and step-bounded WFA2 for HiFi/ONT long reads. |
| **Threats** | • **Skepticism of AI Rewrites:** Community warnings (Nils Homer, Clay McLeod, Kai Blin, `polars-bio`) against "convincingly correct" ports make strict differential CI gating non-negotiable.<br/>• **Deep-Learning & Haplotype Callers:** `DeepSomatic` (*Nat. Biotechnol.* 2025), `ClairS`/`ClairS-TO` (*Nat. Commun.* 2025), `Octopus`, `varlociraptor`, and `sawfish` raise precision expectations in homopolymers and FFPE.<br/>• **Upstream `VarDictJava` EOL (`0886261`):** With `VarDictJava` officially archived (`#396`), `VarDictRust` must steward post-1.0 specifications transparently via [`ROADMAP.md`](../ROADMAP.md). |

---

### 8.2 Comparative Architecture: Fulcrum Genomics, PacBio & Greenfield Rust Callers

| Engineering Dimension | Fulcrum Genomics (`fgumi` v0.9, `perbase` v1.4, `ferro-hgvs` v1.0) | PacBio & Greenfield Rust Suite (`sawfish` v2.3, `trgt` v5.1, `HiPhase` v1.7, `aardvark` v1.1, `varlociraptor` v8.12, `kanpig` v2.0, `modkit` v0.6, `oarfish` v0.11) | `VarDictRust` (`v1.0.7`) Status & Gap Analysis |
| :--- | :--- | :--- | :--- |
| **HTS I/O & Record Representation** | `fgumi` uses `noodles 0.115` + `fgumi-raw-bam` (`RawRecord` zero-copy byte slices) + `noodles-bgzf` (`libdeflate`) + `flate2` (`zlib-rs`); `perbase` uses `rust-htslib` with `thread_local!` `IndexedReader` and `LruCache` reference contigs. | `aardvark`, `oarfish`, and `piscem` use **HTSlib-free `noodles`** (`oarfish` inspects raw packed 4-byte CIGAR words on undecoded `noodles::bam::Record` buffers). `fibertools-rs` had to emergency-pin `hts-sys = "=2.2.0"` due to `bindgen` ABI breakage in `bam1_core_t`. | **Strong Foundation (`main` / `v1.0.8`):** HTSlib-free `noodles 0.118` I/O + native `noodles-bgzf 0.52` (`libdeflate` C backend + `zlib-rs`) avoids `hts-sys` bindgen bugs, with 32-shard `BGZF_BLOCK_CACHE`, LUT sequence decoding, and zero-copy read slices (`PR #98`, `PR #104`). **Remaining Gap:** Raw-record pre-filtering before `BamRecord` materialization (`v2.0`). |
| **Sequence & Alignment Kernels** | `fgumi-dna` packs 4 bases/byte (2-bit DNA, 32 bases/`u64`) for `XOR` + `popcnt` Hamming distance; `fgumi-simd-fastq` uses runtime-dispatched AVX2 / NEON. | `sawfish`, `trgt`, `HiPhase`, and `aardvark` use 2-piece affine **Wavefront Alignment (`WFA2` / `DWFALite`)** + POA (`spoa_rs`). Crucially, `sawfish` caps WFA at `CONSENSUS_SCORE_MAX_WFA_STEPS = 60` (falling back to ~30 µs DP) with a 6 GB global alignment semaphore (`CONTIG_ALIGNMENT_MEMORY_BUDGET_BYTES`). `varlociraptor` uses Myers bit-parallel edit distance to shrink `PairHMM` to an 8 bp band (`EDIT_BAND = 4`). | **Partial:** `find_35_match` uses 64-bit SWAR fast-reject and `realign.rs` reuses scratch buffers with zero-intermediate flanking builders (`PR #98`–`#107`), while 2-bit DNA `XOR`+`popcnt` and step-bounded `block-aligner` / WFA2 remain slated for `v2.0`. |
| **Concurrency, Backpressure & Disk Prefetching** | `fgumi-pipeline-core` uses a 7-stage priority work-stealing scheduler with `sysinfo` cgroup memory limits (`--max-memory auto`) and a 10s deadlock watchdog; `perbase` overlaps superchunk compute and bounded `crossbeam` channel draining via `rayon::join`. | `sawfish` splits threads $\lceil 0.3 N \rceil$ I/O scan vs $0.7 N$ speculative breakpoint assembly, and issues `libc::posix_fadvise(POSIX_FADV_WILLNEED)` on merged BGZF chunks `8 * threads` ahead of `lowest_unfinished`. `modkit` uses a `return_mem` channel so worker threads absorb `drop()` deallocation latency. | **Strong (`main` / `v1.0.8`):** `runner.rs` bounds in-flight tasks (`4 * threads` + `Condvar`), `SpillableVcfAccumulator` spills above 8 MiB, and `scan_multi_tile_max_read_lengths` now runs in parallel across `rayon` (`PR #102`); `DEFAULT_UNTILED_SPAN_MAX = 1_000_000` still leaves 100 kbp–1 Mbp intervals untiled by default. |
| **VCF Output, Indexing & File Safety** | `fgumi-bam-io` performs checked `close(2)` syscalls (preventing `File::drop` from swallowing NFS/Lustre writeback errors); `perbase` handles `ErrorKind::BrokenPipe` cleanly on stdout. | `aardvark` writes `.vcf.gz` via `noodles::bgzf::io::MultithreadedWriter` (`threads.clamp(1, 4)`) and generates `.tbi` purely in Rust (`noodles::vcf::fs::index` + `noodles::tabix::fs::write`). `trgt` reuses a single cleared VCF record buffer across loci. | **Complete (`main` / `v1.0.8`):** Handles `BrokenPipe` cleanly on both `stdout` and `stderr`, writes multithreaded BGZF `.vcf.gz` via `--output` + pure-Rust `--tbi` / `--csi` indexes, and stages regular files via `StagedFileTarget` with `File::sync_all` + input-collision guards (`PR #100`–`#105`). |
| **Testing, Telemetry & Parity Governance** | `fgumi` uses `proptest` (with committed regression seeds), `rstest`, `cargo nextest`, `miri`, and built-in `--pipeline-stats`. `Bismark` v3.1 freezes default v1.0 byte parity while gating v2+ features behind **never-silent** opt-in flags + concordance gates. | `sawfish` (`run.stats.json`) and `HiPhase` (`--summary-file`) emit structured JSON telemetry and high-depth excluded BEDs (`max.sv.depth.bed`); `paraphase` codifies legacy Python parity and intentional low-depth divergences in unit tests. | **Strong Parity & Telemetry (`main` / `v1.0.8`):** 8-gate differential harness + `UPSTREAM.toml` (`VarDictJava`, `VarDict`, and advisory `vardictcpp 3`) + `--stats-json` telemetry + `--parity` / `--compat vardictjava-1.8.3` pinning (`PR #100`–`#107`); `proptest` interval-conservation suites remain a candidate expansion. |

---

### 8.3 Clinical Oncology Literature & Community Pain-Point Synthesis

1. **Where VarDict Wins Clinically (Lai et al., *NAR* 2016; SEQC2, *Nat. Biotechnol.* 2021; Guille et al., *Brief. Bioinform.* 2025):**
   * **Complex Indels & MNVs in Actionable Cancer Genes:** Rescues soft-clipped reads (`CigarModifier` + `VariationRealigner` + `StructuralVariantsProcessor`) to call composite indels and 20–250 bp frameshifts in *BRCA1/2*, *EGFR* (exon 19 deletions / exon 20 insertions in **AURA/FLAURA/ADAURA** ctDNA trials, *J. Thorac. Oncol.* 2024), *TP53*, and *FLT3*-ITD that haplotype/column callers split or miss.
   * **Ultra-Deep Low-VAF Sensitivity Without Downsampling:** Evaluates every read at $1,000\times\text{–}50,000\times$ depth down to $0.1\%$ VAF (Sandmann et al., *Sci. Rep.* 2017).
   * **Native 8-Column BED Amplicon Calling (`-a`):** Masks primer binding sites and flags single-amplicon PCR artifacts (`AMPBIAS`) without external trimming tools.
   * **Simultaneous Somatic, Germline & Candidate LOH Classification:** Heuristic `STATUS` classification in `postprocess/somatic.rs` (`StrongSomatic`, `LikelySomatic`, `Germline`, `StrongLOH`, `LikelyLOH`, `Deletion`), paired with Fisher exact strand-bias and somatic p-value annotations (`SBF`, `ODDRATIO`, `SSF`), identifies somatic variants alongside germline and candidate LOH events in a single pass.
2. **Documented Scientific & Community Pain Points (Bian et al., 2018; Chen et al., 2020; Zhang et al., *DeepFilter* 2023; GitHub / Biostars / Reddit):**
   * **Low-VAF Precision / FDR Without Post-Filtering:** Raw calling at `-f 0.01` yields elevated false positives in homopolymers (`MSI >= 5`) and repeats (on 30× **germline WGS** `HG002`, raw $F_1 = 0.6622$ improves to **$F_1 = 0.9699$** at `QUAL >= 128`; low-VAF somatic and ctDNA panels require separate depth- and VAF-aware calibration).
   * **Missing `F1R2`/`F2R1` Orientation Bias & UMI Consensus:** Lacks native FFPE deamination (`C>T`/`G>A`) / 8-oxoG (`C>A`/`G>T`) read-pair orientation modeling (vs. Mutect2 `LearnReadOrientationModel` or Pisces `VQR`) and native duplex UMI collapsing (vs. `fgumi`).
   * **Upstream Output Quirks Replicated in `v1.0.x`:**
     * `VD > DP` and `AF > 1` (`VarDictJava#92`, `#242`, `#316`, `#323`, `#386`, `#400`) when soft-clip realignment adds supporting reads (`VD`) without incrementing locus `DP`.
     * Single-winner multi-allelic collapsing (`UB-7`, `VarDict#35`, `VarDictJava#65`) compounded by the `-A` flag collision (`sv_min_af` in `VarDictJava` vs `$opt_A` all-alleles in `var2vcf_valid.pl`).
     * Missing automatic 3/5/6-column BED detection: default 13-column UCSC `refGene` settings (`-c 3 -S 7 -E 8 -s 10 -e 11 -g 13`, with 4-column standard BED and 8-column amplicon BED auto-detected) fail with `VardictError::Parse` on 3-, 5-, and 6-column UCSC BED files (`regions.rs:641-652`) when `-c 1 -S 2 -E 3 -z` is omitted (`VarDictJava#49`, `#54`).
     * NovaSeq 2-color dark-cycle poly-G mate-pileup stalls (`sv-lazy-mate-pileup`, `VarDictJava#370`).

---

### 8.4 Prioritized Future Work Matrix (`Tier 1` $\to$ `Tier 2` $\to$ `Tier 3`) & Tracked GitHub Epics (`#90`–`#96`)

To eliminate end-user `--fix` flag sprawl while permanently preserving `VarDictJava v1.8.3` reproducibility, post-1.0 evolution uses a **Two-Phase `--parity` vs `--modern` SemVer Model** ([`#94`](https://github.com/palakpsheth/VarDictRust/issues/94)):
- **`v1.0.x` (Current Baseline):** Parity-only by default; all default CLI invocations are **100% `VarDictJava v1.8.3` byte-identical**, and explicit **`--parity`** / **`--compat vardictjava-1.8.3`** landed on `main` ahead of `v1.0.8` (`PR #102`) so pipelines can forward-compatibly pin v1.8.3 parity today (`--modern` remains rejected until `v1.1.0`).
- **`v1.1.x`:** Default CLI invocations remain **100% `VarDictJava v1.8.3` byte-identical** (`--parity` implied by default). Introduces **`--modern`** (opting into the curated bundle of clinical/spec fixes [`#95`](https://github.com/palakpsheth/VarDictRust/issues/95) and parity-breaking speedups [`#34`](https://github.com/palakpsheth/VarDictRust/issues/34) `PB-1..PB-5`) alongside **`--parity`** (`--compat vardictjava-1.8.3`).
- **`v2.0.0+`:** **`--modern`** becomes the default out-of-the-box behavior, while **`--parity`** reproduces the rigorously tested `VarDictJava v1.8.3` byte-for-byte oracle forever with a single flag.

| Tier & Target Release | Priority Item & Tracked Epic | Subsystem / Location | Summary of Change & Verification Gate |
| :--- | :--- | :--- | :--- |
| **Tier 1 (`v1.0.x`)**<br/>*Zero-Output-Diff Engine, Perf & UX* | **1.1 `fastpath-modify-cigar`** ([`#90`](https://github.com/palakpsheth/VarDictRust/issues/90) — **landed on `main` ahead of `v1.0.8`, `PR #98`, `PR #102`, `PR #104` & `PR #105`**) | `cigar/modify.rs:408-806` | On single-`Match` (`150M`, >85% of WGS reads) and pure `S`/`M` reads, apply endpoint mismatch/clipping checks directly in-place on `record.cigar` with zero-copy `seq`/`qual` slices and `modify_cigar_preserving_buffers` scratch reuse, bypassing the inner 25-regex indel/clip cascade and ASCII string round-trip when unclipped (~11.5% CPU savings). Gate: endpoint-mismatch parity fixtures + `pr_check.sh` 0 diffs. |
| **Tier 1 (`v1.0.x`)** | **1.2 `parallel-max-rl-scan` & `tile-span-100k`** ([`#90`](https://github.com/palakpsheth/VarDictRust/issues/90) — **parallel pre-scan landed on `main` ahead of `v1.0.8`, `PR #102`**) | `lib.rs:69, 412-893, 1677, 2368` | Parallelize per-tile `max_read_length` pre-scan across `rayon` workers with prefix-max fold (re-running only tiles impacted by an upstream `max_read_length` step increase; `DEFAULT_UNTILED_SPAN_MAX = 1_000_000` retained until `tile-span-100k` parity tuning). |
| **Tier 1 (`v1.0.x`)** | **1.3 `direct-loaded-region` & `zero-alloc-realign`** ([`#90`](https://github.com/palakpsheth/VarDictRust/issues/90) — **landed on `main` ahead of `v1.0.8`, `PR #98`–`#107`**) | `reference.rs`, `cigar/parse.rs`, `realign.rs`, `sv.rs`, `variations/mod.rs` | Share `Arc<[u8]>` `Reference` spans (`RegionSpanSlot`), inline `Reference::get_base` primary-region lookup, and pack `[u64; 3]` `SeedKey` tables; inline single-allele `VariationEntries` + `JavaHashMap` new-key delta log, Somatic `RefSiteMode`/`CompactRefSite` elision, `O(W_mate log W_mate)` SV `partial_pipeline` delta log, `join_ref_for_5_lgins` / `join_ref_for_3_lgins` zero-intermediate builders, `i32::MAX` overflow saturation, and `FxHashMap` non-iterated lookup maps. |
| **Tier 1 (`v1.0.x`)** | **1.4 `bgzf-cache-scale` & `tuple-dedup-key`** ([`#90`](https://github.com/palakpsheth/VarDictRust/issues/90) — **landed on `main` ahead of `v1.0.8`, `PR #98`, `PR #104`, `PR #105` & `PR #107`**) | `bam.rs` | Expand `BGZF_BLOCK_CACHE` ($8 \to 32$ shards / pooled buffers) to eliminate 32-thread `fetch_block Arc` eviction churn; add LUT/vectorized BAM sequence decode + zero-copy `seq`/`qual` read slices + checked CIGAR span overflow guards; replace `-t` `HashSet<String>` with `FxHashSet<DupKey>`; construct `pub(crate) build_noodles_region` directly. |
| **Tier 1 (`v1.0.x`)** | **1.5 `audit-fisher-stream` & `audit-somatic-dedup`** ([`#90`](https://github.com/palakpsheth/VarDictRust/issues/90), [`#91`](https://github.com/palakpsheth/VarDictRust/issues/91) — **landed on `main` ahead of `v1.0.8`, `PR #98`–`#107`**) | `stats/fisher.rs`, `stats/f80.rs`, `output/vcf.rs`, `lib.rs` | Eliminate redundant `support: Vec<i32>` and heap `Vec<f64>` allocations in `FisherExact` via algebraic `support_len()`, 4,096-entry LRU cache, `f80_sum_iter`, `build_single_vcf_event` `--pass-only` Fisher/uniroot short-circuit in `output/vcf.rs`, and `StackBuf` / `write_perl_num_str` formatting (`#90`); apply `SplitVarDeduper` (`200 kbp` cross-tile `.var` ownership window, `PR #99`), single-buffer `SingleVcfRecord` / `PairedVcfRecord` consolidation, `SpillableVcfAccumulator::insert` `IndexMap::get_mut(chr)` fast-path, deterministic `-i` splice sorting, and streaming per-tile `BamFeed` pileup (`#91`). |
| **Tier 1 (`v1.0.x`)** | **1.6 CLI & File-Safety Hardening** ([`#93`](https://github.com/palakpsheth/VarDictRust/issues/93) — **landed on `main` ahead of `v1.0.8`, `PR #98`, `PR #100`, `PR #102`, `PR #104`–`#107`**) | `cli.rs`, `main.rs`, `regions.rs`, `runner.rs`, `Cargo.toml` | Sync `cli.rs` `MIT AND GPL-2.0-or-later` license string and MSRV `1.91.0` + `crates.io` metadata; warn on 5/6-column BEDs lacking cols 7/8; apply `StagedFileTarget` atomic staging + `File::sync_all` + input-collision guards across `--output`/`--tbi`/`--csi`/`--stats-json`; support `-T`/`--trim`, `-n`, `-VS SILENT\|LENIENT`, clustered `-pz`/`-uk`/`-ua`, `RAYON_NUM_THREADS`, GRCh38 `HLA-*` colon-contig `-R`/BED parsing (`regions::parse_region_flag`), and all 17 consolidated `var2vcf` CLI flags covering all 16 live Perl post-filter/threshold options plus `--completions` and `--manpage` (`89` CLI flags across `7` functional groups; see [`docs/migration.md`](migration.md)); auto-set 64 MiB worker stack on `.cram` inputs; add `UPSTREAM.toml`. |
| **Tier 2 (`v1.1.0`)**<br/>*Cleanroom MIT, `--parity`/`--modern` & UX* | **2.1 Cleanroom `stats/rmath/` $\to$ Pure `MIT`** ([`#92`](https://github.com/palakpsheth/VarDictRust/issues/92)) | `stats/rmath/`, `stats/fisher.rs:343-478`, `Cargo.toml:9` | Replace the 5 GPL-2.0 R `nmath` files (~350 LoC) with a cleanroom Loader (2000) / Apache-2.0 saddle-point implementation verified to $\varepsilon \le 10^{-6}$, restoring pure `MIT` licensing across Cargo and Bioconda. Because R 4.6's 5-decimal rounding (`1e-5` step vs `1e-6` gate) flips at dyadic midpoints (`35/259` Commons Math 3 and `20/259` exact-rational flips vs `0/259` `stats/rmath/` on `r_fisher_oracle.tsv`), cleanroom MIT Fisher math ships under `v1.1.0` `--modern` while `--parity` preserves exact R 4.6 equivalence. |
| **Tier 2 (`v1.1.0`)** | **2.2 Native `--output out.vcf.gz` + `.tbi`/`.csi`, Consolidated `var2vcf` Flags & `--stats-json`** ([`#93`](https://github.com/palakpsheth/VarDictRust/issues/93) — **landed early on `main` ahead of `v1.0.8`, `PR #100`, `PR #102`, `PR #104`–`#107`**) | `crates/vardictrust/src/cli.rs`, `main.rs` | Direct multithreaded BGZF `.vcf.gz` output + pure-Rust `.tbi` (`--tbi` / `--index`, with actionable `--csi` hint above `512 Mbp`) and `.csi` (`--csi`) index generation (`END=`-aware for structural variants) via `noodles::bgzf::io::MultithreadedWriter` + `noodles::tabix` / `noodles::csi` (`aardvark` pattern), all 17 consolidated `var2vcf` CLI flags covering all 16 live Perl post-filter/threshold options (`--no-end`, `--pass-only`, `--paired-strict`, `--vcf-pval`, `--min-splitreads`, `--no-pstd-filter`, `--vcf-af`, `--vcf-min-depth`, `--vcf-min-vd`, `--vcf-qmean`, `--vcf-mapq`, `--vcf-pmean`, `--vcf-sn`, `--vcf-max-nm`, `--vcf-gt-freq`, `--vcf-msi`, `--vcf-cluster-dist`), `-y / --verbose` progress, and `--stats-json` run telemetry (`sawfish`/`fgumi` pattern). |
| **Tier 2 (`v1.1.0`)** | **2.3 `--parity` / `--modern` Mode & Clinical Fixes** ([`#94`](https://github.com/palakpsheth/VarDictRust/issues/94), [`#95`](https://github.com/palakpsheth/VarDictRust/issues/95)) | [`ROADMAP.md`](../ROADMAP.md) | `--parity` (`--compat vardictjava-1.8.3`) landed ahead of `v1.0.8` (`PR #102`); `v1.1.0` adds `--modern` (`#94`, bundling `#34` `PB-1..PB-5` + `#95` clinical fixes: `vd-dp-reconcile`, `--all-alleles` / `ub7`, `ub5`, `eqx-read-length`, `sv-lazy-mate-pileup`, `378` paired amplicon, `402`, `342`, `292`, `369`, `var2vcf-chrx-drop`) with never-silent VCF header + `stderr` notices (`Bismark` v3.1 pattern). |
| **Tier 2 (`v1.1.0`)** | **2.4 Clinical Presets (`--preset`)** ([`#94`](https://github.com/palakpsheth/VarDictRust/issues/94)) | `crates/vardictrust/src/cli.rs` | Add `--preset wgs\|wes\|amplicon\|ctdna` bundling 3-col BED coordinates (`-c 1 -S 2 -E 3 -z`), soft `FILTER` tags (`#62`), germline WGS `QUAL >= 128` calibration (`HG002` 30× germline $F_1 = 0.9699$), and separate depth/VAF-aware somatic & dilution calibration for `ctdna` / tumor-normal panels. |
| **Tier 3 (`v2.0`)**<br/>*Next-Gen Clinical & Rust-Native* | **3.1 Streaming UMI/Duplex & `F1R2`/`F2R1` Orientation Bias** ([`#96`](https://github.com/palakpsheth/VarDictRust/issues/96)) | `bam.rs`, `cigar/parse.rs`, `postprocess/` | Default to `--modern` (`--parity` preserves v1.8.3 forever); tile-local coordinate-sorted UMI & duplex consensus calling (`--umi-tag RX --duplex`, `fgumi` pattern) + Beta-Binomial `F1R2`/`F2R1` read-orientation bias model (`--orientation-bias-model`) for FFPE & liquid biopsy ctDNA. |
| **Tier 3 (`v2.0`)** | **3.2 2-Bit DNA, Step-Bounded WFA2 & Long-Read Mode** ([`#96`](https://github.com/palakpsheth/VarDictRust/issues/96)) | `realign.rs`, `sv.rs`, `cigar/parse.rs` | 2-bit DNA `XOR`+`popcnt` matching (`fgumi-dna`), Myers band-shrinking (`varlociraptor`), step-bounded 2-piece affine WFA2 (`sawfish` 60-step cap + 6 GB global alignment semaphore), rolling ring-buffer pileup (`[PosSlot; 16384]`), and `--preset hifi\|ont` (`SA:Z` split reads, >125 bp SV clips, `MM`/`ML` & `HP`/`PS` tag preservation). |

