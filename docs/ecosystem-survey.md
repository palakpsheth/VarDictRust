# Ecosystem Survey — AI-assisted Rust rewrites of bioinformatics tools

**Author / Maintainer:** Palak Sheth <palakpsheth@gmail.com>  
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

| Dimension | VarDictRust v1.0.0 | VarDictJava v1.8.3 | GATK Mutect2 (v4.5+) | Illumina Strelka2 (v2.9+) | Octopus (v0.7+) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Language / Runtime** | Single static Rust binary (`musl`/`glibc`) | Java 8 + R 3.x + Perl 5 pipeline | Java JVM (GATK framework) | C++11 + Python wrapper | C++14 native |
| **Memory Footprint** | Bounded ($O(1)$ halo tiling, $\le 2.9\text{ GB}$) | Unbounded ($O(N)$ depth, OOM prone) | High (JVM heap 8–16 GB) | Moderate (4–8 GB) | High (Haplotype graph 8–16 GB) |
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
| **Thread & partition invariance** | Formally verified in [`thread_determinism.rs`](../crates/vardictrust/tests/thread_determinism.rs). |
| **Upstream bugs, manually verified** | Documented in [`docs/upstream-bugs.md`](upstream-bugs.md) ([UB-1](upstream-bugs.md#ub-1--p-spurious-all-zero-record-at-position-00) through [UB-7](upstream-bugs.md#ub-7-multi-allelic-site-single-winner-collapsing-var2vcf_validpl135-205)). |
| **Credit & AI transparency** | Comprehensive attribution and disclosures in [`README.md`](../README.md) and [`benchmarks/evidence/`](../benchmarks/evidence/). |

---

## 5. Architectural Hardening & Verification Matrix

| Area | Historical Gap | Engine Hardening Measure | Verification Evidence |
| :--- | :--- | :--- | :--- |
| **Hostile CIGAR Fuzzing** | Crash risks in CIGAR modifier | Bounds-checked CIGAR manipulation in `crates/vardictrust-core/src/cigar/` | `crates/vardictrust/tests/robustness.rs` covering PR #402 corrupt topology |
| **CLI Flag Parity** | Incomplete CLI flag documentation | Explicit 60-flag Diátaxis reference categorized into 6 functional groups | Section 4 in [`README.md`](../README.md) |
| **Release Artifact Packaging** | Ungated release binaries | PIE static binary validation gate comparing musl to glibc outputs | [`harness/bin_equiv.sh`](../harness/bin_equiv.sh), [`musl-ci.json`](../benchmarks/evidence/musl-ci.json) |
| **Biological Truth Concordance** | Lack of physical truth-set evaluation | Full multi-cohort RTG `vcfeval` truth concordance suite | 5 cohorts in [`docs/validation/omnibus-empirical-validation.md`](validation/omnibus-empirical-validation.md) |
| **Upstream Drift Defense** | Unpinned upstream dependency risk | Pinned reference oracle to VarDictJava v1.8.3 (`commit 0886261`) | Sealed environment specification in [`harness/setup_oracle.sh`](../harness/setup_oracle.sh) |

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
