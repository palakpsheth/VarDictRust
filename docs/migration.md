# Migrating from VarDictJava + R + Perl to VarDictRust

`vardictrust` consolidates the three-stage [`VarDictJava` v1.8.3](https://github.com/AstraZeneca-NGS/VarDictJava/tree/0886261) + [`VarDict` R/Perl](https://github.com/AstraZeneca-NGS/VarDict/tree/009e017) pipeline (`VarDictJava` $\to$ `teststrandbias.R` / `testsomatic.R` $\to$ `var2vcf_valid.pl` / `var2vcf_paired.pl`) into a **single statically-compiled binary** with zero JVM, R, or Perl runtime dependencies.

On valid inputs, `vardictrust` produces VCF 4.2 (`VCFv4.2`) and intermediate `.var` output **byte-for-byte identical** (within floating-point tolerance $\varepsilon \le 10^{-6}$ after R 5-decimal banker's rounding) to `VarDictJava` v1.8.3 (`0886261`) + `VarDict` (`009e017`). See [`DIVERGENCES.md`](../DIVERGENCES.md) and [`docs/validation/omnibus-empirical-validation.md`](validation/omnibus-empirical-validation.md) for empirical parity receipts.

---

## 1. Side-by-Side Pipeline Translation

### Single-Sample / Germline Calling

**Legacy 3-process pipe (`VarDictJava` + R + Perl on a standard 0-based UCSC BED):**

```bash
VarDict \
  -G ref.fa \
  -f 0.01 \
  -N SAMPLE \
  -b sample.bam \
  -c 1 -S 2 -E 3 -g 4 -z 1 \
  targets.bed \
  | teststrandbias.R \
  | var2vcf_valid.pl -N SAMPLE -E -f 0.01 \
  > sample.vcf
```

**Single-binary `vardictrust` equivalent (stdout VCF or direct BGZF + Tabix):**

```bash
# Uncompressed VCF to stdout (drop-in replacement):
vardictrust \
  -G ref.fa \
  -f 0.01 \
  -N SAMPLE \
  -b sample.bam \
  -c 1 -S 2 -E 3 -g 4 -z \
  --no-end \
  targets.bed \
  > sample.vcf

# Or write multithreaded BGZF + Tabix (.tbi) directly with zero samtools/bgzip/tabix pipes:
vardictrust \
  -G ref.fa \
  -f 0.01 \
  -N SAMPLE \
  -b sample.bam \
  -c 1 -S 2 -E 3 -g 4 -z \
  --no-end \
  --output sample.vcf.gz \
  --tbi \
  targets.bed
```

### Paired Tumor/Normal Somatic Calling

**Legacy 3-process pipe (`VarDictJava` + `testsomatic.R` + `var2vcf_paired.pl` on a standard 0-based UCSC BED):**

```bash
VarDict \
  -G ref.fa \
  -f 0.01 \
  -N "TUMOR|NORMAL" \
  -b "tumor.bam|normal.bam" \
  -c 1 -S 2 -E 3 -g 4 -z 1 \
  targets.bed \
  | testsomatic.R \
  | var2vcf_paired.pl -N "TUMOR|NORMAL" -M -f 0.01 \
  > somatic.vcf
```

**Single-binary `vardictrust` equivalent:**

```bash
vardictrust \
  -G ref.fa \
  -f 0.01 \
  -N "TUMOR|NORMAL" \
  -b "tumor.bam|normal.bam" \
  -c 1 -S 2 -E 3 -g 4 -z \
  --paired-strict \
  --output somatic.vcf.gz \
  --tbi \
  targets.bed
```

### Targeted 8-Column Amplicon Calling (`-a` / `AMPBIAS`)

```bash
# 8-column amplicon BED (chr, gene_start, gene_end, gene, score, strand, insert_start, insert_end)
vardictrust \
  -G ref.fa \
  -f 0.01 \
  -N SAMPLE \
  -b amplicon.bam \
  -a 10:0.95 \
  amplicon_8col.bed \
  > amplicon.vcf
```

### Intermediate `.var` Tabular Output (`--emit-intermediate`)

When migrating custom downstream scripts that consume raw `.var` tables instead of VCF:

- **`vardictrust ... --emit-intermediate`**: emits the unannotated Java `.var` schema (`36` columns in Simple mode, `38` in Amplicon mode, `55` in Somatic mode — identical to running `VarDictJava` alone without R).
- **`vardictrust ... --emit-intermediate --fisher`**: runs the in-engine R-exact Fisher strand-bias / somatic tests and emits the R-annotated `.var` schema (`38` columns in Simple mode, `40` in Amplicon mode, `61` in Somatic mode — identical to `VarDictJava | teststrandbias.R` or `VarDictJava | testsomatic.R`).

---

## 2. `var2vcf_valid.pl` & `var2vcf_paired.pl` Flag Rosetta Stone

Because `vardictrust` fuses [`VarDictJava`](https://github.com/AstraZeneca-NGS/VarDictJava/blob/0886261/src/main/java/com/astrazeneca/vardict/CmdParser.java) and the Perl [`var2vcf_valid.pl`](https://github.com/AstraZeneca-NGS/VarDict/blob/009e017/var2vcf_valid.pl) / [`var2vcf_paired.pl`](https://github.com/AstraZeneca-NGS/VarDict/blob/009e017/var2vcf_paired.pl) post-filters into one CLI, Perl single-letter flags that collided with Java flags (`-S`, `-E`, `-M`, `-P`, `-T`, `-d`, `-v`, `-q`, `-Q`, `-p`, `-o`, `-m`, `-F`, `-I`, `-c`) map to explicit `--vcf-*` / long flags in Group 2 (`VCF Post-Filtering (var2vcf)`):

| Legacy Perl Flag (`var2vcf_valid.pl` / `var2vcf_paired.pl`) | `vardictrust` Flag | Single-Sample Default (`valid`) | Paired Somatic Default (`paired`) | Notes & Exact Parity Semantics |
| :--- | :--- | :--- | :--- | :--- |
| `-N <NAME>` (`TUMOR\|NORMAL` in paired) | `-N <NAME>` / `--name <NAME>` | BAM stem | `<bam1>\|<bam2>` stems | Shared directly between pileup engine and VCF `#CHROM` sample columns. |
| `-E` (`var2vcf_valid.pl` only) | `--no-end` | Emit `END=` | Never emits `END=` | Omits `END=<INT>` in single-sample VCF `INFO` (paired mode already omits `END=`). |
| `-S` | `--pass-only` | Emit all (`FILTER` tagged) | Emit all (`FILTER` tagged) | Emits only `PASS` variants; non-`PASS` sites short-circuit before expensive VCF Fisher formatting when safe. |
| `-M` (`var2vcf_paired.pl` only) | `--paired-strict` | N/A | Off (germline rescue on) | Disables germline/LOH rescue and restricts paired output to strict somatic candidates (`var2vcf_paired.pl -M`). |
| `-P <FLOAT>` (`var2vcf_paired.pl`) | `--vcf-pval <FLOAT>` | N/A | `0.05` | Maximum somatic Fisher p-value (`SSF`) threshold and `P<FLOAT>` `FILTER` tag in paired VCF output. |
| `-P 0` (`var2vcf_valid.pl`) | `--no-pstd-filter` | Filter `PSTD=0` on | N/A | Disables the single-position read filter (`PSTD == 0`) in single-sample VCF output. |
| `-T <INT>` (`var2vcf_valid.pl`) | `--min-splitreads <INT>` | `1` | N/A | Minimum `SPLITREAD` support required for structural variants in single-sample VCF output. |
| `-f <FLOAT>` | `--vcf-af <FLOAT>` (or `-f <FLOAT>`) | Inherits `-f` (`0.02` if `-f 0`) | Inherits `-f` | Overrides VCF `FILTER`/`GT` minimum AF independently of pileup `-f`. Passing `--vcf-af 0` explicitly sets `0.0` (`defined()` semantics) in both modes; when `--vcf-af` is omitted, single-sample `var2vcf_valid.pl` preserves Perl's falsy fallback (`-f 0` $\to$ `0.02`; pass `-f 0.0` or `--vcf-af 0` for zero single-sample VCF AF cutoff), whereas paired somatic `var2vcf_paired.pl` uses `defined($opt_f)` so `-f 0` stays `0`. |
| `-d <INT>` | `--vcf-min-depth <INT>` | `3` | `5` | Minimum total read depth (`DP`) before `d<INT>` soft-filter. |
| `-v <INT>` | `--vcf-min-vd <INT>` | `2` | `3` | Minimum variant read depth (`VD`) before `v<INT>` soft-filter. |
| `-q <FLOAT>` | `--vcf-qmean <FLOAT>` | `22.5` | `22.5` | Minimum mean base quality (`QUAL`) before `q<FLOAT>` soft-filter. |
| `-Q <FLOAT>` | `--vcf-mapq <FLOAT>` | `10` | `0` | Minimum mean mapping quality (`MQ`) before `Q<FLOAT>` soft-filter. |
| `-p <FLOAT>` | `--vcf-pmean <FLOAT>` | `8` | `8` | Minimum mean read position (`PMEAN`) before `p<FLOAT>` soft-filter. |
| `-o <FLOAT>` | `--vcf-sn <FLOAT>` | `1.5` | `1.5` | Minimum signal-to-noise ratio (`SN`) before `SN<FLOAT>` soft-filter. |
| `-m <FLOAT>` | `--vcf-max-nm <FLOAT>` | `5.25` | `5.25` | Maximum mean mismatches (`NM`) before `NM<FLOAT>` soft-filter. |
| `-F <FLOAT>` | `--vcf-gt-freq <FLOAT>` | `0.2` | `0.2` | Homozygous `1/1` genotype margin (`1 - AF < FLOAT`). |
| `-I <FLOAT>` | `--vcf-msi <FLOAT>` | `12` | `12` | Maximum non-monomer microsatellite repeat count (`MSI`) before `MSI<FLOAT>` filter. |
| `-c <INT>` | `--vcf-cluster-dist <INT>` | `0` (disabled) | `0` (disabled) | Adjacent variant cluster window in bp (`Cluster<INT>bp` filter across single-sample, amplicon, and paired somatic streams). |
| `-A` (`var2vcf_valid.pl` all-alleles) | *Deferred to `v1.1.0` `--all-alleles`* | Single-winner (`UB-7`) | Single-winner (`UB-7`) | **Important:** `-A <INT>` on `vardictrust` is `VarDictJava`'s `sv_min_af` (`[default: 4]`), **not** `var2vcf_valid.pl -A`! Multi-allelic emission (`UB-7`) is tracked for `v1.1.0` (`--modern` / `--all-alleles`; see [`ROADMAP.md`](../ROADMAP.md)). |

---

## 3. Four Common Migration Gotchas

### 1. `-h` Prints the `.var` Column Header — Use `-H`, `-?`, or `--help` for Help

Following [`CmdParser.java`](https://github.com/AstraZeneca-NGS/VarDictJava/blob/0886261/src/main/java/com/astrazeneca/vardict/CmdParser.java#L67), `-h` / `--header` sets `printHeader = true` (emitting the `.var` column header row in `--emit-intermediate` mode). To print the 7-group CLI help screen, run:

```bash
vardictrust -H
# or:
vardictrust --help
```

### 2. BED Column Indices & 0-Based Coordinate Indexing (`-c 1 -S 2 -E 3 -z`)

Following [`RegionBuilder.java`](https://github.com/AstraZeneca-NGS/VarDictJava/blob/0886261/src/main/java/com/astrazeneca/vardict/RegionBuilder.java#L48-L68), `VarDictJava`'s default non-amplicon table layout is 13-column UCSC `refGene` (`-c 3 -S 7 -E 8 -s 10 -e 11 -g 13`, **1-based inclusive**; 8-column amplicon BEDs are auto-detected separately):

- **4-column UCSC BED (`chr`, `start_0based`, `end_1based`, `gene`) invoked *without* `-c`:** both `VarDictJava` and `vardictrust` auto-detect `-c 1 -S 2 -E 3 -g 4` **and** 0-based half-open coordinates (`-z 1`) when `start <= end`.
- **3-, 5-, and 6-column 0-based UCSC BEDs — and 4-column 0-based BEDs whenever `-c 1` is passed explicitly:** you **must** pass **`-c 1 -S 2 -E 3 -z`** (or `-z 1`, plus `-g 4` when column 4 is the region/gene name)! Passing `-c 1` without `-z` disables 4-column auto-detection and leaves `zero_based = false` (**1-based inclusive**), silently shifting every region start coordinate by 1 bp; omitting `-c 1 -S 2 -E 3` on `< 11`-column non-4-column BEDs fails fast with `VardictError::Parse` expecting at least 11 `refGene` columns. Note that `vardictrust` accepts both bare `-z` and `-z 1` even immediately before the positional `BED` path, whereas legacy `VarDictJava` requires explicit `-z 1` there (Commons CLI `hasOptionalArgs(1)` otherwise consumes the `BED` path as `-z`'s numeric value). Omit `-z` only when your region file uses 1-based inclusive coordinates.

### 3. Single-Dash Legacy Java Flags vs GNU `--` Flags

- **Accepted with one or two dashes (`VarDictJava` Commons-CLI compatibility):** `-th` / `--threads`, `-UN` / `--UN`, `-fisher` / `--fisher`, `-chimeric` / `--chimeric`, `-deldupvar` / `--deldupvar`, `-nosv` / `--nosv`, `-splice` / `--splice`, `-header` / `--header`, `-dedup` / `--dedup`, `-debug` / `--debug`, `-verbose` / `--verbose`, `-amplicon` / `--amplicon`, `-DP` / `--DP` / `--default-printer`, `-adaptor` / `--adaptor`, `-mfreq` / `--mfreq`, `-nmfreq` / `--nmfreq`, `-insert-size` / `--insert-size`, `-insert-std` / `--insert-std`, `-ref-extension` / `--ref-extension`, `-crispr` / `--crispr`, `-trim` / `--trim`, and `-VS` / `--VS`.
- **Require double-dash `--` (`vardictrust`-native flags):** `--output`, `--tbi` (`--index`), `--csi`, `--stats-json`, `--parity`, `--compat`, `--emit-intermediate`, `--tile-size`, `--completions`, `--manpage`, `--no-end`, `--pass-only`, `--paired-strict`, `--min-splitreads`, `--no-pstd-filter`, and `--vcf-*`. Passing a single-dash typo (for example `-pass-only` or `-completions`) fails fast with `Unknown single-dash option '-pass-only'; did you mean '--pass-only'?` rather than being misinterpreted as a short-flag cluster.

### 4. Default Thread Count (`-th` / `RAYON_NUM_THREADS`)

`VarDictJava` runs single-threaded unless `-th <INT>` is passed. Because `vardictrust` guarantees worker-count-invariant deterministic output and bounded $O(1)$ per-tile memory, omitting `-th` defaults to `RAYON_NUM_THREADS` (if set) or all available logical cores. Inside shared HPC/SLURM/Nextflow allocations, pin threads explicitly via `-th <INT>` (for example `-th 8`).

---

## 4. Direct Compressed I/O, Indexing, Telemetry & Shell Ergonomics

`vardictrust` eliminates external wrapper plumbing around compressed inputs, BGZF outputs, and container streams:

- **CRAM & Companion Index Discovery:** Pass `-b sample.cram -G ref.fa` directly. Both standard (`sample.bam.bai`, `sample.cram.crai`, `ref.fa.fai`) and HTSJDK companion (`sample.bai`, `sample.crai`, `ref.fai`, plus symlink targets) index layouts are discovered automatically.
- **BGZF-Compressed Reference FASTA & BED (`v1.0.8+`):** Pass `-G ref.fa.gz` when indexed with `bgzip` + `samtools faidx ref.fa.gz` (`.fai` + `.gzi` present; plain non-BGZF `gzip` FASTAs fail fast with a `bgzip` + `samtools faidx` remediation hint), and pass `.bed.gz`, `.bed.bgz`, or process-substitution FIFOs (`<(zcat targets.bed.gz)`) directly as the positional `BED` input.
- **Direct BGZF + Tabix/CSI Output & Stream Safety:** Use `--output variants.vcf.gz --tbi` (or `--csi` when contigs exceed `512 Mbp`) for atomic staged file output with `fsync` verification, or `--output -` / `/dev/stdout` in streaming container pipelines.
- **Run Telemetry (`--stats-json <PATH>`):** Emit structured JSON metrics (`mode`, `threads`, `regions`, `tiles`, `elapsed_ms`, `peak_rss_kb`) alongside variant calling.
- **Forward-Compatible Parity Pinning (`--parity`):** Pass `--parity` (or `--compat vardictjava-1.8.3`) in production workflows to lock 100% `VarDictJava` v1.8.3 parity across both `v1.x` and future `v2.0+` releases.
- **Shell Completions & Manual Page (`v1.0.8+`):**
  ```bash
  # Generate shell completions (installed automatically in Bioconda environments):
  vardictrust --completions bash > ~/.local/share/bash-completion/completions/vardictrust
  vardictrust --completions zsh  > ~/.zfunc/_vardictrust
  vardictrust --completions fish > ~/.config/fish/completions/vardictrust.fish

  # View the standalone roff manual page:
  vardictrust --manpage | man -l -
  ```

---

## 5. Unsupported Legacy `VarDictJava` Options (Fail-Loud Policy)

Per [`AGENTS.md §12`](../AGENTS.md#12-unsupported-features-error-loudly), `vardictrust` **never** silently ignores unimplemented legacy flags:

| Legacy Flag | Upstream Behavior | `vardictrust` v1.0 Behavior & Recommended Alternative |
| :--- | :--- | :--- |
| `-C` | Strip `chr` prefix / treat chromosome names as bare integers (deprecated upstream). | Rejected loudly (`-C not implemented in v1.0; use VarDictJava`). Match BED/region contig names directly to your BAM/FASTA header (`chr1` vs `1`). |
| `-Z <FLOAT>` / `--downsample <FLOAT>` | Nondeterministic `Math.random()` read subsampling (`CmdParser.java:125-127`). | Rejected loudly (`-Z not implemented in v1.0; use VarDictJava`). Halo tiling bounds memory without downsampling; if subsampling is required for benchmarking, run `samtools view -s <SEED.FRAC>` beforehand for deterministic reproducibility. |
| `-VS STRICT` | HTSJDK strict SAM validation aborting on benign BAM tag quirks. | Rejected loudly (`-VS STRICT not implemented in v1.0; use VarDictJava`). `-VS SILENT` and `-VS LENIENT` are accepted and match `noodles` lenient decoding. |
