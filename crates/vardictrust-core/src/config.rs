//! Engine, var2vcf post-filter, and [`CompatProfile`] compatibility configuration.

/// Opt-in VCF post-filter numeric override carrying an optional normalized CLI token
/// (<= 64 ASCII bytes) so `##FILTER=<ID=...>` and record `FILTER` tags reproduce Perl
/// dualvar interpolation verbatim (e.g. `--vcf-qmean 20.0` -> `q20.0`, programmatic
/// `VcfThreshold::new(20.0)` -> `q20`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfThreshold<T> {
    pub value: T,
    pub raw: Option<String>,
}

impl<T> VcfThreshold<T> {
    #[inline]
    pub const fn new(value: T) -> Self {
        Self { value, raw: None }
    }

    #[inline]
    pub fn with_raw(value: T, raw: String) -> Self {
        Self {
            value,
            raw: Some(raw),
        }
    }
}

impl<T> From<T> for VcfThreshold<T> {
    #[inline]
    fn from(value: T) -> Self {
        Self::new(value)
    }
}

/// Per-sample amplicon concordance counts (`GDAMP`, `TLAMP`, `NCAMP`, `AMPFLAG`)
/// for single-sample and paired somatic amplicon VCF emission (`AstraZeneca-NGS/VarDictJava#378`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AmpliconConcordance {
    pub gdamp: i32,
    pub tlamp: i32,
    pub ncamp: i32,
    pub ampflag: bool,
}

/// Compatibility profile selector modeled on `bwa-mem3`'s `--compat` target enum.
///
/// Default across `v1.x` is [`Self::VardictJava183`] (100% VarDictJava v1.8.3 / `var2vcf` parity).
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum CompatTarget {
    /// Exact VarDictJava v1.8.3 (`0886261`) + `var2vcf_*.pl` (`009e017`) output parity (`--compat=vardictjava`).
    #[default]
    VardictJava183,
    /// Compatibility off (`--compat=off`): enables all 32 [`CompatProfile`] accuracy, scale, SV, amplicon, VCF schema, and `--eqx` fixes and stamps `##VarDictRustCompat=off` in the VCF header.
    Off,
}

impl std::fmt::Display for CompatTarget {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::VardictJava183 => f.write_str("vardictjava"),
            Self::Off => f.write_str("off"),
        }
    }
}

impl CompatTarget {
    /// Return the resolved [`CompatProfile`] for this target.
    #[inline]
    pub const fn profile(self) -> CompatProfile {
        match self {
            Self::VardictJava183 => COMPAT_VARDICTJAVA_1_8_3,
            Self::Off => COMPAT_OFF,
        }
    }
}

/// Resolved compatibility behavior switches (32 independent correctness, scale, SV, amplicon, VCF schema, and `--eqx` flags) derived from [`CompatTarget`].
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct CompatProfile {
    pub target: CompatTarget,
    /// `#116`: Bijective VCF 4.2 percent-encoding on `INFO/SAMPLE=...` (`%`, `;`, `=`, `,`, whitespace/control, `""` -> `-`, `-` -> `%2D`, `.` -> `%2E`).
    pub sanitize_sample_info: bool,
    /// `#342`: Normalize `<dup-N>` (`mlen < 100`) to `<DUP>` with non-negative `i64` span (including `-L <= 100`) and unify homozygous `genotype1` `<DUP>` span (`2 * SVFLANK + dup_val`, `ToVarsBuilder.java:619` vs `:700`).
    pub standard_dup_symbolic_alt: bool,
    /// `#292`: Shift 3'-moved (`-3`) deletion/insertion `REF`/`ALT` alleles to match shifted coordinates.
    pub shift3_3prime_alleles: bool,
    /// `#369`: Restrict emitted variants to merged requested target intervals (`start_pos >= 1`).
    pub strict_region_bounds: bool,
    /// `UB-1`/`UB-2`: Suppress `-p` pileup `<R,R>`/`N/N` zero-depth sentinel rows in Simple and Amplicon modes.
    pub suppress_pileup_sentinels: bool,
    /// Populate `DUPRATE` from live `u64` `RecordPreprocessor` duplicate/total read counts instead of `0`.
    pub live_duprate: bool,
    /// Preserve both unprefixed and `chr`-prefixed sex/mito contigs (`X`+`chrX`, `Y`+`chrY`, `MT`+`chrM`) in VCF output.
    pub preserve_mixed_sex_contigs: bool,
    /// `PR #402`: Absorb negative trailing `M` (`rm < 0`) into `rdoff_corr` during homopolymer multi-indel merges in `CigarModifier` (returning `Err(JavaThrows)` immediately if `rdoff_corr <= 0`).
    pub cigar_homopolymer_merge: bool,
    /// `eqx-read-length`: Treat `CigarKind::SequenceMatch` (`=`) and `CigarKind::SequenceMismatch` (`X`) equivalently to `CigarKind::Match` (`M`) across read-length, coalescing, realignment, dedup, and prescan paths.
    pub eqx_cigar_read_length: bool,
    /// `UB-5`: Preserve single-sample symbolic SVs (`<DUP>`, `<INV>`, `<DEL>`) with `spanpairs > 0` when `splitreads == 0` unless `--min-splitreads` was explicitly overridden.
    pub preserve_discordant_svs: bool,
    /// `UB-6`: Pass pre-decrement cluster end `mark_pe` to `mark_dup_sv` in `svfdup` when `tbp != 0 && tbp < pe`.
    pub mark_dup_sv_inclusive_end: bool,
    /// `sv-lazy-mate-pileup`: Gate `find_inv_sub_paths` mate-window loads on local soft-clip presence (`>= 12 bp`) and defer `partial_pipeline_paths` until `find_match_rev` + gap validation succeed (`sv.rs` `find_inv_sub_paths` / `StructuralVariantsProcessor.java:627-628`).
    pub sv_lazy_mate_pileup: bool,
    /// `rounding-underflow-tie`: Bypass OpenJDK `DecimalFormat` `DigitList` underflow truncation (`±5e-{4,5,6,9}`) so exact binary IEEE-754 `f64` values round to `±10^-d`.
    pub exact_binary_half_even: bool,
    /// #92 / #386 / #400 / VarDict #113: Enforce VD <= DP, AF <= 1.0, ADJAF <= 1.0, and HICNT <= HICOV (HIAF <= 1.0) across all variant types.
    pub reconcile_vd_dp: bool,
    /// UB-3 (#282) / #34 PB-1: Iterate variant and post-processing position maps in ascending genomic coordinate order instead of JavaHashMap bucket order.
    pub coordinate_sorted_vars: bool,
    /// #34 PB-2 (VariationRealigner.java:2006-2042, :2016-2024): Omit Java's out-of-bounds while-loop sentinel +1 mismatch and +1 length increment when rem5 < rem3 in find_35_match.
    pub exact_35_overlap: bool,
    /// #34 PB-4 (VariationRealigner.java:2056-2065): Match exact deletion CIGAR operator length (len == dlen) in no_passing_reads instead of Java's getCigarString().contains(dlen + "D") suffix match.
    pub exact_cigar_dlen_match: bool,
    /// `joachimwolff/VarDictJava#4`, `clintval/vartovcf#139`, `clintval/vartovcf#149` (`SomaticPostProcessModule.java:156-226`): Continue evaluating subsequent candidate alleles at a multi-allelic locus when an earlier candidate fails `is_good_var`.
    pub somatic_continue_after_bad_var: bool,
    /// `AstraZeneca-NGS/VarDictJava#240`, `AstraZeneca-NGS/VarDictJava#241`, `AstraZeneca-NGS/VarDict#119` (`VariationUtils.java:136-146`): Pass directional `5`/`3` soft-clip orientation across all consensus assembly sites and thread `AdaptorSet` through `adjust_mnp` and `adj_snv`.
    pub adaptor_all_consensus_dirs: bool,
    /// `AstraZeneca-NGS/VarDictJava#398`, `AstraZeneca-NGS/VarDictJava#225`, `AstraZeneca-NGS/VarDictJava#229`, `AstraZeneca-NGS/VarDict#111`, `clintval/vartovcf#107`, `clintval/vartovcf#192` (VCFv4.2 §1.4.5/§5.4): Emit canonical 1-bp-upstream-anchored symbolic `<DEL>`, `<DUP>`, and `<INV>` alleles with signed `SVLEN` and paired somatic `END`/`SVTYPE`/`SVLEN` `INFO` tags.
    pub canonical_sv_vcf_anchoring: bool,
    /// `AstraZeneca-NGS/VarDictJava#378`, `AstraZeneca-NGS/VarDictJava#155`, `AstraZeneca-NGS/VarDict#109`, `MHH-Bioinformatics-Hematology/vardictcpp#1`, `palakpsheth/VarDictRust#95`: Enable bounded-memory paired tumor/normal somatic amplicon calling (`-b "tumor|normal" -a`) in both `.var` and VCF modes.
    pub paired_somatic_amplicon: bool,
    /// `AstraZeneca-NGS/VarDictJava#341`, `AstraZeneca-NGS/VarDictJava#103`: Permanently strip outer positive-length hard-clip/pad (`H`/`P`) CIGAR operators without advancing alignment start coordinates, and treat interior `P` as non-consuming.
    pub cigar_strip_outer_hard_clips: bool,
    /// `bcbio/bcbio-nextgen#1704`: Floor subtracted mean read position, base quality, mapping quality, mismatch count (`NM`), and VCF `QUAL` at `0` (`>= 0.0`).
    pub nonnegative_subtracted_means: bool,
    /// `AstraZeneca-NGS/VarDictJava#316`, `clintval/vartovcf#120`, `clintval/vartovcf#142`: Reset all auxiliary `Variation` count/sum accumulators alongside `cnt = 0` and reconcile `fwd + rev == VD` and `0 <= hicnt <= VD` (`locnt = VD - hicnt`).
    pub reconcile_strand_and_hq_sums: bool,
    /// `AstraZeneca-NGS/VarDict#168`, `AstraZeneca-NGS/VarDictJava#323`: Enforce `DP >= RD + VD` (`total_pos_cov >= ref_fwd + ref_rev + pos_coverage`) across all variant and somatic placeholder builders.
    pub reconcile_rd_vd_dp: bool,
    /// `AstraZeneca-NGS/VarDictJava#267`, `AstraZeneca-NGS/VarDictJava#287`, `AstraZeneca-NGS/VarDict#137`, `bcbio/bcbio-nextgen#3321`: Reduce IUPAC ambiguity codes in `ALT` alleles outside symbolic `<...>` tags and skip non-pileup `REF == ALT` records.
    pub sanitize_iupac_alt_alleles: bool,
    /// `AstraZeneca-NGS/VarDictJava#328`, `AstraZeneca-NGS/VarDict#157`: Respect `-L` (`--min-sv-len`) by materializing bounded explicit sequence alleles for split-refined SVs shorter than `-L` and filtering discordant-only clusters shorter than `-L`.
    pub sv_respect_min_len: bool,
    /// `AstraZeneca-NGS/VarDictJava#330`: Hydrate and emit the reference pileup row at `-p` positions when all surviving variants shift away via deletions or `adj_complex`.
    pub pileup_preserve_shifted_ref: bool,
    /// `clintval/vartovcf#109`, `clintval/vartovcf#147`: Recompute `Variant::var_type` after `adj_complex` trims a complex allele into an `SNV`, `Insertion`, or `Deletion`.
    pub recompute_vartype_after_adj_complex: bool,
    /// `clintval/vartovcf#101`, `AstraZeneca-NGS/VarDict#75`, `AstraZeneca-NGS/VarDict#76`, `EBIvariation/vcf-validator#99`: Emit single-sample `FORMAT/AD` with `Number=R` cardinality (`RD` when `ALT == "."`, `RD,VD` including `RD,0` when `ALT != "."`).
    pub strict_vcf_ad_number_r: bool,
    /// `clintval/vartovcf#105`: Omit `INFO/END` on non-symbolic VCF records by default (retaining `INFO/END` strictly on symbolic `<...>` SV alleles unless `--no-end` is set).
    pub omit_nonsv_vcf_end: bool,
    /// `clintval/vartovcf#99`, `clintval/vartovcf#54`, `AstraZeneca-NGS/VarDictJava#344`, `AstraZeneca-NGS/VarDict#186`: Canonicalize unphased `1/0` to `0/1` and promote concrete non-dot `0/0` carrier genotypes with `VD > 0` to `0/1`.
    pub canonical_vcf_gt: bool,
}

/// Default `v1.x` VarDictJava v1.8.3 / `var2vcf_*.pl` parity profile (`--compat=vardictjava`, all 32 [`CompatProfile`] fixes disabled).
pub const COMPAT_VARDICTJAVA_1_8_3: CompatProfile = CompatProfile {
    target: CompatTarget::VardictJava183,
    sanitize_sample_info: false,
    standard_dup_symbolic_alt: false,
    shift3_3prime_alleles: false,
    strict_region_bounds: false,
    suppress_pileup_sentinels: false,
    live_duprate: false,
    preserve_mixed_sex_contigs: false,
    cigar_homopolymer_merge: false,
    eqx_cigar_read_length: false,
    preserve_discordant_svs: false,
    mark_dup_sv_inclusive_end: false,
    sv_lazy_mate_pileup: false,
    exact_binary_half_even: false,
    reconcile_vd_dp: false,
    coordinate_sorted_vars: false,
    exact_35_overlap: false,
    exact_cigar_dlen_match: false,
    somatic_continue_after_bad_var: false,
    adaptor_all_consensus_dirs: false,
    canonical_sv_vcf_anchoring: false,
    paired_somatic_amplicon: false,
    cigar_strip_outer_hard_clips: false,
    nonnegative_subtracted_means: false,
    reconcile_strand_and_hq_sums: false,
    reconcile_rd_vd_dp: false,
    sanitize_iupac_alt_alleles: false,
    sv_respect_min_len: false,
    pileup_preserve_shifted_ref: false,
    recompute_vartype_after_adj_complex: false,
    strict_vcf_ad_number_r: false,
    omit_nonsv_vcf_end: false,
    canonical_vcf_gt: false,
};

/// Compatibility-off profile (`--compat=off`) enabling all 32 [`CompatProfile`] accuracy, scale, symbolic SV (`<DUP>`, `<INV>`, `<DEL>`), paired somatic amplicon, VCF schema, and `--eqx` fixes.
pub const COMPAT_OFF: CompatProfile = CompatProfile {
    target: CompatTarget::Off,
    sanitize_sample_info: true,
    standard_dup_symbolic_alt: true,
    shift3_3prime_alleles: true,
    strict_region_bounds: true,
    suppress_pileup_sentinels: true,
    live_duprate: true,
    preserve_mixed_sex_contigs: true,
    cigar_homopolymer_merge: true,
    eqx_cigar_read_length: true,
    preserve_discordant_svs: true,
    mark_dup_sv_inclusive_end: true,
    sv_lazy_mate_pileup: true,
    exact_binary_half_even: true,
    reconcile_vd_dp: true,
    coordinate_sorted_vars: true,
    exact_35_overlap: true,
    exact_cigar_dlen_match: true,
    somatic_continue_after_bad_var: true,
    adaptor_all_consensus_dirs: true,
    canonical_sv_vcf_anchoring: true,
    paired_somatic_amplicon: true,
    cigar_strip_outer_hard_clips: true,
    nonnegative_subtracted_means: true,
    reconcile_strand_and_hq_sums: true,
    reconcile_rd_vd_dp: true,
    sanitize_iupac_alt_alleles: true,
    sv_respect_min_len: true,
    pileup_preserve_shifted_ref: true,
    recompute_vartype_after_adj_complex: true,
    strict_vcf_ad_number_r: true,
    omit_nonsv_vcf_end: true,
    canonical_vcf_gt: true,
};

impl Default for CompatProfile {
    #[inline]
    fn default() -> Self {
        COMPAT_VARDICTJAVA_1_8_3
    }
}

/// Parse a `--compat <TARGET>` string into [`CompatTarget`].
///
/// # Errors
///
/// Returns [`crate::types::VardictError::Unsupported`] (`Err(VardictError::Unsupported(...))`)
/// when `raw` is not a recognized compatibility target (`vardictjava`, `vardictjava-1.8.3`,
/// `1.8.3`, or `off`).
pub fn parse_compat_target(raw: &str) -> Result<CompatTarget, crate::types::VardictError> {
    let trimmed = raw.trim();
    if trimmed.eq_ignore_ascii_case("vardictjava")
        || trimmed.eq_ignore_ascii_case("vardictjava-1.8.3")
        || trimmed == "1.8.3"
    {
        Ok(CompatTarget::VardictJava183)
    } else if trimmed.eq_ignore_ascii_case("off") {
        Ok(CompatTarget::Off)
    } else {
        Err(crate::types::VardictError::Unsupported(format!(
            "--compat {raw} not supported; expected 'vardictjava' (or 'vardictjava-1.8.3' / '1.8.3') or 'off'"
        )))
    }
}

/// Central engine and `var2vcf` post-filter configuration struct.
///
/// Preserves all VarDictJava v1.8.3 (`Configuration.java`) and `var2vcf_*.pl` defaults;
/// constructed programmatically via [`Config::default()`] or from CLI arguments via
/// `Cli::to_config()` in the `vardictrust` binary crate.
#[derive(Debug, Clone)]
pub struct Config {
    pub exception_counter: std::sync::Arc<std::sync::atomic::AtomicUsize>,
    /// Reference FASTA path.
    pub reference: String,
    /// BAM/CRAM path(s). Somatic: "tumor.bam|normal.bam".
    pub bam: String,
    /// Allele frequency threshold (default: 0.01).
    pub freq: f64,
    /// Sample name(s).
    pub sample: String,
    /// True when `-N` was explicitly passed on the CLI (`Configuration.sample` != null before `initResources`).
    pub sample_explicit: bool,
    /// Region string (if -R provided).
    pub region: Option<String>,
    /// BED file path (positional).
    pub bed: Option<String>,
    /// Delimiter for BED files (default: "\t")
    pub delimiter: String,
    /// BED column indices (0-based internally).
    pub col_chr: usize,
    pub col_start: usize,
    pub col_end: usize,
    pub col_seg_start: usize,
    pub col_seg_end: usize,
    pub col_gene: Option<usize>,
    /// Track if column for chromosome was explicitly set
    pub column_for_chromosome: Option<usize>,
    /// Thread count.
    pub threads: usize,
    /// Local realignment (0 or 1).
    pub realign: u8,
    /// Minimum variant reads (default: 2).
    pub min_reads: u32,
    /// Base quality threshold (default: 22.5).
    pub base_qual: f64,
    /// Mapping quality filter.
    pub map_qual: Option<u32>,
    /// Read position filter (default: 5.0).
    pub read_pos: f64,
    /// Max mismatches (default: 8).
    pub max_mismatches: u32,
    /// Minimum match+insertion read length; reads shorter are ignored (Java -M, default: 0).
    pub minmatch: u32,
    /// Low frequency threshold for somatic STATUS determination (default: 0.05, -V flag).
    pub lofreq: f64,
    /// Disable SV calling.
    pub no_sv: bool,
    /// SV minimum length (default: 1000).
    pub sv_min_len: u32,
    /// Emit intermediate .var instead of VCF.
    pub emit_intermediate: bool,
    /// Include Fisher p-value/OR columns.
    pub fisher: bool,
    /// Tile size for bounded memory (None = auto).
    pub tile_size: Option<u32>,
    /// Amplicon mode string (e.g. "10:0.95").
    pub amplicon: Option<String>,
    /// Single-primer amplicon mode (`-sp` / `--amplicon-single-primer`, ports `AstraZeneca-NGS/VarDictJava#303`).
    pub amplicon_single_primer: bool,

    // -- Additional flags added in T027-T032 --
    /// Extension for BED region (default: 0).
    pub extension: i32,
    /// Number of bases to extend for reads (default: 2).
    pub extension_nuc: i32,
    /// Shift indels to 3'-most alternative alignment (-3 / moveIndelsTo3).
    pub three_prime: bool,
    /// Minimum # of reads to call with (default: 2).
    pub min_reads_bias: u32,
    /// qratio threshold (default: 1.5).
    pub qratio: f64,
    /// MapQ for output (`Configuration.outputMapq: Double`, default: 0.0).
    pub out_mapq: f64,
    /// Indel size (default: 50).
    pub indel_size: u32,
    /// Zero-based BED coordinates.
    pub zero_based: bool,
    /// True when `-z` was explicitly supplied on the command line (`Configuration.isZeroBasedDefined()`).
    pub zero_based_defined: bool,
    /// SAM flag filter (default: 0x504).
    pub flag_filter: u32,
    /// Raw `-F` CLI string (`Configuration.samfilter: String`, default: `"0x504"`).
    pub samfilter_raw: String,
    /// SV parameters.
    pub sv_read_len: u32,
    pub sv_flanking: u32,
    pub sv_min_af: u32,
    pub max_region_len: u32,
    /// Mono/di frequency thresholds.
    pub mfreq: f64,
    pub nmfreq: f64,
    /// Chimeric read handling.
    pub chimeric: bool,
    /// Delete/dup variant handling.
    pub deldupvar: bool,
    /// Splice mode.
    pub splice: bool,
    /// Pileup-only mode.
    pub pileup: bool,
    /// Debug mode (-D): append per-read debug columns to intermediate output.
    pub debug: bool,
    /// Debug print target.
    pub debug_print: Option<String>,
    /// CRISPR parameters.
    pub crispr_cutting_site: Option<i32>,
    pub crispr_flanking_size: Option<i32>,
    /// Column header output.
    pub column_header: bool,
    /// Dedup flag.
    pub dedup: bool,
    /// List of adaptor sequences.
    pub adaptor: Vec<String>,
    /// Verbose mode.
    pub verbose: bool,
    /// The literal -f value for the VCF writer's $Freq (var2vcf_valid.pl /
    /// var2vcf_paired.pl get it from their own command line, independent of the
    /// engine's pileup override that sets `freq` to -1, CmdParser.java:134-136).
    pub af_threshold: f64,
    /// True when `-f` was given as a non-`"0"` zero string (e.g. `-f 0.0` or `-f 0.00`),
    /// which is truthy in Perl (`var2vcf_valid.pl:21` `my $Freq = $opt_f ? $opt_f : 0.02;`)
    /// whereas `-f 0` (`"0"`) is falsy and falls back to `0.02`.
    pub af_threshold_truthy_zero: bool,
    /// Normalized raw CLI token supplied to `-f` / `--af-threshold` for verbatim `f$Freq`
    /// VCF header/FILTER formatting (`None` when `-f` was not explicitly passed;
    /// `AstraZeneca-NGS/VarDict@009e017`, `MHH-Bioinformatics-Hematology/vardictcpp#4`).
    pub af_threshold_raw: Option<String>,
    /// unique mode alignment enabled (-u)
    pub unique_mode_alignment_enabled: bool,
    /// unique mode second in pair enabled (-UN)
    pub unique_mode_second_in_pair_enabled: bool,
    /// Omit END tag in VCF output (matches var2vcf_valid.pl -E)
    pub omit_end_tag: bool,
    /// True only when `--no-end` (or `--noend`) was explicitly passed on the CLI
    /// (`false` in `Config::default()`, so symbolic SVs under `canonical_sv_vcf_anchoring`
    /// emit canonical VCF `INFO/END` unless `--no-end` / `--noend` was explicitly requested).
    pub vcf_no_end_explicit: bool,
    /// Include Ns in depth calculation (-K)
    pub include_n: bool,
    /// Opt-in output file path (`--output <PATH>`; plain file or `.gz`/`.bgz` BGZF).
    pub output_path: Option<String>,
    /// Write a Tabix (.tbi) index alongside a BGZF-compressed --output .vcf.gz file (--tbi / --index).
    pub index_output: bool,
    /// Opt-in JSON execution telemetry path (`--stats-json <PATH>`).
    pub stats_json_path: Option<String>,
    /// Trim bases after N bases from 5' end (`-T` / `--trim-bases-after`, default: 0 = disabled).
    pub trim_bases_after: i32,
    /// Sample name regular expression (`-n` / `--sample-regexp`, with outer `/` stripped).
    pub sample_regexp: Option<String>,
    /// Write a CSI (`.csi`) index alongside a BGZF-compressed `--output .vcf.gz` file (`--csi`).
    pub index_csi: bool,
    /// Emit only PASS variants in VCF output (`--pass-only`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -S`).
    pub vcf_pass_only: bool,
    /// Disable germline rescue in paired somatic VCF filtering (`--paired-strict`, ports `var2vcf_paired.pl -M`).
    pub vcf_paired_strict: bool,
    /// Maximum p-value threshold for high-confidence somatic calls (`--vcf-pval`, default: 0.05, ports `var2vcf_paired.pl -P`).
    pub vcf_pval: f64,
    /// Normalized raw CLI token supplied to `--vcf-pval` for verbatim `P$PVAL` formatting (`None` when unset).
    pub vcf_pval_raw: Option<String>,
    /// Minimum split reads required for structural variants in single-sample VCF (`--min-splitreads`, default: 1, ports `var2vcf_valid.pl -T`).
    pub vcf_min_splitreads: i32,
    /// True when `--min-splitreads` was explicitly passed on the CLI (`var2vcf_valid.pl -T`).
    pub vcf_min_splitreads_explicit: bool,
    /// Enable the PSTD=0 single-position read filter in single-sample VCF output (`true` by default; `--no-pstd-filter` sets `false`, ports `var2vcf_valid.pl -P 0`).
    pub vcf_pstd_filter: bool,
    /// Opt-in VCF allele-frequency filter/GT threshold override (`--vcf-af`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -f`).
    pub vcf_af: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum total read depth override (`--vcf-min-depth`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -d`).
    pub vcf_min_depth: Option<VcfThreshold<i32>>,
    /// Opt-in VCF minimum variant depth override (`--vcf-min-vd`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -v`).
    pub vcf_min_vd: Option<VcfThreshold<i32>>,
    /// Opt-in VCF minimum mean base quality override (`--vcf-qmean`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -q`).
    pub vcf_qmean: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum mean mapping quality override (`--vcf-mapq`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -Q`).
    pub vcf_mapq: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum mean read-position override (`--vcf-pmean`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -p`).
    pub vcf_pmean: Option<VcfThreshold<f64>>,
    /// Opt-in VCF minimum signal-to-noise (`qratio`) override (`--vcf-sn`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -o`).
    pub vcf_sn: Option<VcfThreshold<f64>>,
    /// Opt-in VCF maximum mean mismatches (`NM`) override (`--vcf-max-nm`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -m`).
    pub vcf_max_nm: Option<VcfThreshold<f64>>,
    /// Opt-in VCF homozygous `1/1` allele-frequency margin override (`--vcf-gt-freq`, `1 - AF < F`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -F`).
    pub vcf_gt_freq: Option<f64>,
    /// Opt-in VCF non-monomer MSI repeat-count threshold override (`--vcf-msi`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -I`).
    pub vcf_msi: Option<VcfThreshold<f64>>,
    /// Opt-in VCF SNV cluster window in bp (`--vcf-cluster-dist`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -c`).
    pub vcf_cluster_dist: Option<VcfThreshold<i32>>,
    /// Compatibility profile (`--compat=vardictjava|off`; default: [`COMPAT_VARDICTJAVA_1_8_3`]).
    pub compat: CompatProfile,
    /// Emit all alleles passing filters at the same locus in VCF output (`--vcf-all-alleles` / `--all-alleles`, ports `var2vcf_valid.pl` / `var2vcf_paired.pl -A`).
    pub vcf_all_alleles: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            exception_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            reference: String::new(),
            bam: String::new(),
            freq: 0.01,
            af_threshold: 0.01,
            af_threshold_truthy_zero: false,
            af_threshold_raw: None,
            sample: "SAMPLE".to_string(),
            sample_explicit: false,
            region: None,
            bed: None,
            delimiter: "\t".to_string(),
            col_chr: 2,
            col_start: 6,
            col_end: 7,
            col_seg_start: 9,
            col_seg_end: 10,
            col_gene: Some(12),
            column_for_chromosome: None,
            threads: num_cpus(),
            realign: 1,
            min_reads: 2,
            base_qual: 22.5,
            map_qual: None,
            read_pos: 5.0,
            max_mismatches: 8,
            minmatch: 0,
            lofreq: 0.05,
            no_sv: false,
            sv_min_len: 1000,
            emit_intermediate: false,
            fisher: false,
            tile_size: None,
            amplicon: None,
            amplicon_single_primer: false,
            extension: 0,
            extension_nuc: 2,
            three_prime: false,
            min_reads_bias: 2,
            qratio: 1.5,
            out_mapq: 0.0,
            indel_size: 50,
            zero_based: false,
            zero_based_defined: false,
            flag_filter: 0x504,
            samfilter_raw: "0x504".to_string(),
            sv_read_len: 300,
            sv_flanking: 100,
            sv_min_af: 4,
            max_region_len: 1200,
            mfreq: 0.25,
            nmfreq: 0.1,
            chimeric: false,
            deldupvar: false,
            splice: false,
            pileup: false,
            debug: false,
            debug_print: None,
            crispr_cutting_site: None,
            crispr_flanking_size: None,
            column_header: false,
            dedup: false,
            adaptor: Vec::new(),
            verbose: false,
            unique_mode_alignment_enabled: false,
            unique_mode_second_in_pair_enabled: false,
            omit_end_tag: true,
            vcf_no_end_explicit: false,
            include_n: false,
            output_path: None,
            index_output: false,
            stats_json_path: None,
            trim_bases_after: 0,
            sample_regexp: None,
            index_csi: false,
            vcf_pass_only: false,
            vcf_paired_strict: false,
            vcf_pval: 0.05,
            vcf_pval_raw: None,
            vcf_min_splitreads: 1,
            vcf_min_splitreads_explicit: false,
            vcf_pstd_filter: true,
            vcf_af: None,
            vcf_min_depth: None,
            vcf_min_vd: None,
            vcf_qmean: None,
            vcf_mapq: None,
            vcf_pmean: None,
            vcf_sn: None,
            vcf_max_nm: None,
            vcf_gt_freq: None,
            vcf_msi: None,
            vcf_cluster_dist: None,
            compat: COMPAT_VARDICTJAVA_1_8_3,
            vcf_all_alleles: false,
        }
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn validate_opt_vcf_f64_threshold(
    slot: &mut Option<VcfThreshold<f64>>,
    field_name: &str,
) -> Result<(), crate::types::VardictError> {
    if let Some(thr) = slot.as_mut() {
        if !thr.value.is_finite() || thr.value < 0.0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid {field_name} '{}': must be a finite number >= 0.0",
                thr.value
            )));
        }
        if thr.value == 0.0 {
            thr.value = 0.0;
        }
        if let Some(ref mut raw) = thr.raw {
            let (norm_val, norm_str) = validate_vcf_nonneg_float_token(raw).map_err(|reason| {
                crate::types::VardictError::Parse(format!(
                    "Invalid {field_name}.raw '{raw}': {reason}"
                ))
            })?;
            if norm_val.to_bits() != thr.value.to_bits() {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid {field_name}.raw '{raw}': parsed value {norm_val} disagrees with {field_name} '{}'",
                    thr.value
                )));
            }
            if raw.as_str() != norm_str {
                *raw = norm_str.to_string();
            }
        }
    }
    Ok(())
}

fn validate_opt_vcf_i32_threshold(
    slot: &mut Option<VcfThreshold<i32>>,
    field_name: &str,
) -> Result<(), crate::types::VardictError> {
    if let Some(thr) = slot.as_mut() {
        if thr.value < 0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid {field_name} '{}': must be an integer >= 0",
                thr.value
            )));
        }
        if let Some(ref mut raw) = thr.raw {
            let (norm_val, norm_str) = validate_vcf_nonneg_i32_token(raw).map_err(|reason| {
                crate::types::VardictError::Parse(format!(
                    "Invalid {field_name}.raw '{raw}': {reason}"
                ))
            })?;
            if norm_val != thr.value {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid {field_name}.raw '{raw}': parsed value {norm_val} disagrees with {field_name} '{}'",
                    thr.value
                )));
            }
            if raw.as_str() != norm_str {
                *raw = norm_str.to_string();
            }
        }
    }
    Ok(())
}

impl Config {
    /// Apply pileup mode (`-p`) parameter overrides matching `Configuration.java:162-165`.
    pub fn apply_pileup_overrides(&mut self) {
        if self.pileup {
            self.freq = -1.0;
            self.min_reads = 0;
        }
    }

    /// Ports `!conf.samfilter.equals("0")` in `SamView.java:101` (`0x100` secondary alignment skip)
    /// and `CigarParser.java:367` (`0x800` supplementary alignment skip).
    /// When `samfilter_raw` was updated from default `"0x504"`, checks `samfilter_raw != "0"`;
    /// otherwise falls back to `flag_filter != 0` so programmatic tests setting `config.flag_filter = 0` work.
    pub fn skip_secondary_and_supplementary(&self) -> bool {
        if self.samfilter_raw != "0x504" {
            self.samfilter_raw != "0"
        } else {
            self.flag_filter != 0
        }
    }

    /// Validates and canonicalizes all VCF post-filter thresholds (`vcf_pval`,
    /// `vcf_min_splitreads`, and all 11 `vcf_*` override fields), canonicalizing `-0.0`/`-0`
    /// and verifying every `raw` CLI token via [`validate_vcf_nonneg_float_token`] /
    /// [`validate_vcf_nonneg_i32_token`] agreement so programmatic `&Config` callers share
    /// identical guards.
    pub fn validate_vcf_thresholds(&mut self) -> Result<(), crate::types::VardictError> {
        if self.af_threshold == 0.0 {
            self.af_threshold = 0.0;
        }
        if let Some(ref mut raw) = self.af_threshold_raw {
            let (f_norm, norm_raw) = validate_vcf_nonneg_float_token(raw).map_err(|reason| {
                crate::types::VardictError::Parse(format!(
                    "Invalid af_threshold_raw '{raw}': {reason}"
                ))
            })?;
            if f_norm.to_bits() != self.af_threshold.to_bits() {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid af_threshold_raw '{raw}': parsed value {f_norm} disagrees with af_threshold '{}'",
                    self.af_threshold
                )));
            }
            if raw.as_str() != norm_raw {
                *raw = norm_raw.to_string();
            }
        }
        if !self.vcf_pval.is_finite() || self.vcf_pval < 0.0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid vcf_pval '{}': must be a finite number >= 0.0",
                self.vcf_pval
            )));
        }
        if self.vcf_pval == 0.0 {
            self.vcf_pval = 0.0;
        }
        if let Some(ref mut raw) = self.vcf_pval_raw {
            let (p_norm, norm_raw) = validate_vcf_nonneg_float_token(raw).map_err(|reason| {
                crate::types::VardictError::Parse(format!("Invalid vcf_pval_raw '{raw}': {reason}"))
            })?;
            if p_norm.to_bits() != self.vcf_pval.to_bits() {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid vcf_pval_raw '{raw}': parsed value {p_norm} disagrees with vcf_pval '{}'",
                    self.vcf_pval
                )));
            }
            if raw.as_str() != norm_raw {
                *raw = norm_raw.to_string();
            }
        }
        if self.vcf_min_splitreads < 0 {
            return Err(crate::types::VardictError::Parse(format!(
                "Invalid vcf_min_splitreads '{}': must be an integer >= 0",
                self.vcf_min_splitreads
            )));
        }

        validate_opt_vcf_f64_threshold(&mut self.vcf_af, "vcf_af")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_qmean, "vcf_qmean")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_mapq, "vcf_mapq")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_pmean, "vcf_pmean")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_sn, "vcf_sn")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_max_nm, "vcf_max_nm")?;
        validate_opt_vcf_f64_threshold(&mut self.vcf_msi, "vcf_msi")?;

        validate_opt_vcf_i32_threshold(&mut self.vcf_min_depth, "vcf_min_depth")?;
        validate_opt_vcf_i32_threshold(&mut self.vcf_min_vd, "vcf_min_vd")?;
        validate_opt_vcf_i32_threshold(&mut self.vcf_cluster_dist, "vcf_cluster_dist")?;

        if let Some(ref mut gt_freq) = self.vcf_gt_freq {
            if !gt_freq.is_finite() || *gt_freq < 0.0 || *gt_freq > 1.0 {
                return Err(crate::types::VardictError::Parse(format!(
                    "Invalid vcf_gt_freq '{gt_freq}': must be a finite number in [0.0, 1.0]"
                )));
            }
            if *gt_freq == 0.0 {
                *gt_freq = 0.0;
            }
        }

        Ok(())
    }
}

/// Ports `java.lang.Integer.decode(String)` (`Configuration.java:155` / `SamView.java:53`):
/// - Rejects empty strings and strings with leading/trailing whitespace.
/// - Handles optional leading `+` or `-` sign.
/// - Detects radix from prefix: `0x` / `0X` / `#` -> 16, leading `0` with remaining length > 1 -> 8, else 10.
/// - Rejects a second sign character after the radix prefix (`0x-1`, `0x+1`, `#-1`).
/// - Parses into `i64` and validates `i32::MIN..=i32::MAX` (so `"0xFFFFFFFF"` is rejected).
pub fn java_integer_decode(s: &str) -> Result<i32, String> {
    if s.is_empty() || s.trim() != s {
        return Err(format!("For input string: \"{s}\""));
    }
    let (negative, rest) = if let Some(r) = s.strip_prefix('-') {
        (true, r)
    } else if let Some(r) = s.strip_prefix('+') {
        (false, r)
    } else {
        (false, s)
    };
    let (radix, digits) =
        if let Some(r) = rest.strip_prefix("0x").or_else(|| rest.strip_prefix("0X")) {
            (16u32, r)
        } else if let Some(r) = rest.strip_prefix('#') {
            (16u32, r)
        } else if let Some(r) = rest.strip_prefix('0').filter(|r| !r.is_empty()) {
            (8u32, r)
        } else {
            (10u32, rest)
        };
    if digits.is_empty() || digits.starts_with('-') || digits.starts_with('+') {
        return Err(format!("For input string: \"{s}\""));
    }
    let mag =
        i64::from_str_radix(digits, radix).map_err(|_| format!("For input string: \"{s}\""))?;
    let val = if negative { -mag } else { mag };
    i32::try_from(val).map_err(|_| format!("For input string: \"{s}\""))
}

/// Ports `org.apache.commons.cli.TypeHandler.createNumber(String).doubleValue()` and
/// `Double.parseDouble` (`CmdParser.java:114`):
/// - If trimmed `s` (after stripping optional trailing `f`/`F`/`d`/`D`) contains `.`, `e`, or `E`:
///   parses as finite `f64` (rejecting `NaN`/`Infinity`).
/// - Otherwise: parses as `Long.valueOf(s)` (strict base-10 signed 64-bit integer, no whitespace)
///   and converts to `f64`.
pub fn parse_java_number_f64(s: &str) -> Result<f64, String> {
    let trimmed = s.trim();
    let core = trimmed
        .strip_suffix(['f', 'F', 'd', 'D'])
        .unwrap_or(trimmed);
    if trimmed.contains('.') || trimmed.contains('e') || trimmed.contains('E') {
        if core.is_empty() {
            return Err(format!("Invalid numeric value '{s}'"));
        }
        let v = core
            .parse::<f64>()
            .map_err(|_| format!("Invalid numeric value '{s}'"))?;
        if !v.is_finite() {
            return Err(format!("Invalid numeric value '{s}'"));
        }
        Ok(v)
    } else {
        if s.is_empty() || s.trim() != s {
            return Err(format!("Invalid numeric value '{s}'"));
        }
        let v = s
            .parse::<i64>()
            .map_err(|_| format!("Invalid numeric value '{s}'"))?;
        Ok(v as f64)
    }
}

/// Ports `((Number) cmd.getParsedOptionValue(...)).intValue()` (`CmdParser.java`).
pub fn parse_java_number_i32(s: &str) -> Result<i32, String> {
    if s.contains('.') {
        let v = parse_java_number_f64(s)?;
        Ok(v as i32)
    } else {
        if s.is_empty() || s.trim() != s {
            return Err(format!("Invalid integer value '{s}'"));
        }
        let v = s
            .parse::<i64>()
            .map_err(|_| format!("Invalid integer value '{s}'"))?;
        Ok(v as i32)
    }
}

/// Unsigned `u32` wrapper around [`parse_java_number_i32`] rejecting negative results.
pub fn parse_java_number_u32(s: &str) -> Result<u32, String> {
    let v = parse_java_number_i32(s)?;
    u32::try_from(v).map_err(|_| format!("Invalid non-negative integer value '{s}'"))
}

/// Unsigned `usize` wrapper around [`parse_java_number_i32`] rejecting negative results.
pub fn parse_java_number_usize(s: &str) -> Result<usize, String> {
    let v = parse_java_number_i32(s)?;
    usize::try_from(v).map_err(|_| format!("Invalid non-negative integer value '{s}'"))
}

/// Unsigned `u8` wrapper around [`parse_java_number_i32`] rejecting out-of-range results.
pub fn parse_java_number_u8(s: &str) -> Result<u8, String> {
    let v = parse_java_number_i32(s)?;
    u8::try_from(v).map_err(|_| format!("Invalid u8 value '{s}'"))
}

fn strip_outer_slashes(raw: &str) -> &str {
    let pat = raw.strip_prefix('/').unwrap_or(raw);
    pat.strip_suffix('/').unwrap_or(pat)
}

const SAMPLE_REGEX_SIZE_LIMIT: usize = 1 << 20;

fn compile_sample_regex(pat: &str) -> Result<regex::Regex, regex::Error> {
    regex::RegexBuilder::new(pat)
        .size_limit(SAMPLE_REGEX_SIZE_LIMIT)
        .build()
}

fn extract_sample_name_clean(bam_path: &str, clean_pat: Option<&str>) -> String {
    static SAMPLE_RE1: std::sync::LazyLock<Option<regex::Regex>> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"([^/\._]+).sorted[^/]*\.bam").ok());
    static SAMPLE_RE2: std::sync::LazyLock<Option<regex::Regex>> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"([^/]+)[_\.][^/]*bam").ok());

    if let Some(pat) = clean_pat {
        if let Ok(re) = compile_sample_regex(pat) {
            if let Some(caps) = re.captures(bam_path) {
                if let Some(m) = caps.get(1) {
                    return m.as_str().to_string();
                }
            }
        }
    } else if let Some(re1) = SAMPLE_RE1.as_ref() {
        if let Some(caps) = re1.captures(bam_path) {
            if let Some(m) = caps.get(1) {
                return m.as_str().to_string();
            }
        }
    }

    if let Some(re2) = SAMPLE_RE2.as_ref() {
        if let Some(caps) = re2.captures(bam_path) {
            if let Some(m) = caps.get(1) {
                return m.as_str().to_string();
            }
        }
    }

    "null".to_string()
}

/// Extract sample name from a BAM file path (or `bam1|bam2` specification),
/// matching `VarDictLauncher.java:75-102` and `Patterns.java:65-67`:
/// 1. `SAMPLE_PATTERN` (`([^/\\._]+).sorted[^/]*\.bam`) or custom `-n` regex (with leading/trailing `/` stripped)
///    applied to the full `bam_path` (`conf.bam.getBamRaw()`).
/// 2. `SAMPLE_PATTERN2` (`([^/]+)[_\.][^/]*bam`) applied to the full `bam_path`.
/// 3. Fallback `"null"` (`VarDictLauncher.java:76`).
pub fn extract_sample_name(bam_path: &str, sample_regexp: Option<&str>) -> String {
    let clean_pat = sample_regexp
        .map(strip_outer_slashes)
        .filter(|s| !s.is_empty());
    extract_sample_name_clean(bam_path, clean_pat)
}

/// Validates a `-n` / `--sample-regexp` CLI pattern after stripping a single
/// leading `'/'` and single trailing `'/'` (`CmdParser.java:88-92`).
///
/// Uses the same linear-time `regex` DFA compiler (`size_limit = 1 MiB`) as
/// `extract_sample_name` and `resolve_sample_names`; non-regular backtracking
/// constructs (lookaround/backreferences) fail fast per `AGENTS.md §12` rather
/// than silently falling back to `"null"`.
pub fn validate_sample_regexp(raw: &str) -> Result<Option<String>, String> {
    let pat = strip_outer_slashes(raw);
    if pat.is_empty() {
        return Ok(None);
    }
    compile_sample_regex(pat).map_err(|e| format!("Invalid sample name regex '{raw}': {e}"))?;
    Ok(Some(pat.to_string()))
}

/// Validates a non-negative VCF threshold `<FLOAT>` token without allocating.
///
/// Enforces `1..=64` trimmed ASCII float bytes (`[0-9.+-eE]`), rejects non-finite numbers
/// (`NaN`/`inf`/overflow) and negative values (including negative non-zero mantissas that
/// underflow `f64` to `-0.0` such as `"-1e-400"`), and canonicalizes signed zeros
/// (`"-0"` -> `(0.0, "0")`, `"-0.0"` / `"+0.0"` -> `(0.0, "0.0")`, `"+0.00"` -> `(0.0, "0.00")`)
/// while preserving verbatim positive formatting (`"1e-5"`, `"0.050"`, `"20.0"`).
pub fn validate_vcf_nonneg_float_token(s: &str) -> Result<(f64, &str), &'static str> {
    if s.is_empty() || s.trim() != s {
        return Err("must be a non-empty trimmed float token");
    }
    if s.len() > 64 {
        return Err("must be <= 64 ASCII bytes");
    }
    if !s
        .bytes()
        .all(|b| matches!(b, b'0'..=b'9' | b'.' | b'+' | b'-' | b'e' | b'E'))
    {
        return Err("must contain only ASCII float characters [0-9.+-eE]");
    }
    let v = s
        .parse::<f64>()
        .map_err(|_| "invalid floating-point syntax")?;
    let neg_nonzero_mantissa = s.starts_with('-')
        && s.split(['e', 'E'])
            .next()
            .is_some_and(|m| m.bytes().any(|b| b.is_ascii_digit() && b != b'0'));
    if !v.is_finite() || v < 0.0 || neg_nonzero_mantissa {
        return Err("must be a finite number >= 0.0");
    }
    if v == 0.0 {
        Ok((0.0, s.strip_prefix(['-', '+']).unwrap_or(s)))
    } else {
        Ok((v, s))
    }
}

/// Inline alias for [`validate_vcf_nonneg_float_token`] retained for `--vcf-pval` callers.
#[inline]
pub fn validate_vcf_pval_token(s: &str) -> Result<(f64, &str), &'static str> {
    validate_vcf_nonneg_float_token(s)
}

/// Validates and normalizes a `--vcf-pval <FLOAT>` (`var2vcf_paired.pl -P`) token into
/// `(f64, String)` for CLI parsing and [`Config`] validation.
pub fn parse_vcf_pval_token(s: &str) -> Result<(f64, String), String> {
    validate_vcf_nonneg_float_token(s)
        .map(|(v, norm)| (v, norm.to_string()))
        .map_err(|reason| format!("Invalid --vcf-pval value '{s}': {reason}"))
}

/// Parses a non-negative VCF threshold `<FLOAT>` flag (`flag`, e.g. `"--vcf-qmean"`) into
/// [`VcfThreshold<f64>`], preserving the normalized raw token for Perl dualvar formatting.
pub fn parse_vcf_threshold_nonneg_f64(s: &str, flag: &str) -> Result<VcfThreshold<f64>, String> {
    let (value, norm) = validate_vcf_nonneg_float_token(s)
        .map_err(|reason| format!("Invalid {flag} value '{s}': {reason}"))?;
    Ok(VcfThreshold::with_raw(value, norm.to_string()))
}

/// Parses a unit-interval `[0.0, 1.0]` VCF float flag (`flag`, e.g. `"--vcf-gt-freq"`).
/// Canonicalizes signed zero (`-0.0` -> `0.0`) and rejects values `> 1.0` or `< 0.0`.
pub fn parse_vcf_unit_f64(s: &str, flag: &str) -> Result<f64, String> {
    let (value, _) = validate_vcf_nonneg_float_token(s)
        .map_err(|reason| format!("Invalid {flag} value '{s}': {reason}"))?;
    if value > 1.0 {
        return Err(format!(
            "Invalid {flag} value '{s}': must be a finite number in [0.0, 1.0]"
        ));
    }
    Ok(value)
}

/// Parses a trimmed base-10 non-negative `i32` (`>= 0`) for `--min-splitreads`.
pub fn parse_vcf_nonneg_i32(s: &str, flag: &str) -> Result<i32, String> {
    if s.is_empty() || s.trim() != s {
        return Err(format!("Invalid {flag} value '{s}'"));
    }
    let v = s
        .parse::<i32>()
        .map_err(|e| format!("Invalid {flag} value '{s}': {e}"))?;
    if v < 0 {
        return Err(format!(
            "Invalid {flag} value '{s}': must be an integer >= 0"
        ));
    }
    Ok(v)
}

/// Validates a non-negative VCF threshold `<INT>` token without allocating.
///
/// Enforces `1..=64` trimmed ASCII integer bytes (`[0-9+-]`), rejects negative values,
/// and canonicalizes `-0` (`"-0"` -> `(0, "0")`) while preserving leading zeros and `+`
/// signs verbatim (`"03"`, `"+02"`, `"010"`, `"00"`, `"+0"`) to match Perl dualvar
/// interpolation in `var2vcf_valid.pl` and `var2vcf_paired.pl`.
pub fn validate_vcf_nonneg_i32_token(s: &str) -> Result<(i32, &str), &'static str> {
    if s.is_empty() || s.trim() != s {
        return Err("must be a non-empty trimmed integer token");
    }
    if s.len() > 64 {
        return Err("must be <= 64 ASCII bytes");
    }
    if !s.bytes().all(|b| matches!(b, b'0'..=b'9' | b'+' | b'-')) {
        return Err("must contain only ASCII integer characters [0-9+-]");
    }
    let v = s.parse::<i32>().map_err(|_| "invalid integer syntax")?;
    if v < 0 {
        return Err("must be an integer >= 0");
    }
    if v == 0 {
        Ok((0, s.strip_prefix('-').unwrap_or(s)))
    } else {
        Ok((v, s))
    }
}

/// Parses a non-negative VCF threshold `<INT>` flag (`flag`, e.g. `"--vcf-min-depth"`) into
/// [`VcfThreshold<i32>`], preserving the normalized raw token for Perl dualvar formatting.
pub fn parse_vcf_threshold_nonneg_i32(s: &str, flag: &str) -> Result<VcfThreshold<i32>, String> {
    let (value, norm) = validate_vcf_nonneg_i32_token(s)
        .map_err(|reason| format!("Invalid {flag} value '{s}': {reason}"))?;
    Ok(VcfThreshold::with_raw(value, norm.to_string()))
}

/// Resolves sample name(s) matching `VarDictLauncher.java:82-86, 206-268` and `var2vcf_paired.pl:45-48`.
/// Splits `bam_raw` ONLY on `'|'` (never on `':'`).
pub fn resolve_sample_names(
    bam_raw: &str,
    explicit_sample: Option<&str>,
    sample_regexp: Option<&str>,
    has_region: bool,
    emit_intermediate_or_splice: bool,
) -> String {
    let clean_pat = sample_regexp
        .map(strip_outer_slashes)
        .filter(|s| !s.is_empty());
    if let Some((bam1, bam2)) = bam_raw.split_once('|') {
        if has_region && emit_intermediate_or_splice {
            let mut base1 = match explicit_sample {
                Some(exp) if clean_pat.is_some() => exp.to_string(),
                Some(exp) => exp.split('|').next().unwrap_or(exp).to_string(),
                None => extract_sample_name_clean(bam_raw, clean_pat),
            };
            if let Some(pat) = clean_pat {
                if let Ok(re) = compile_sample_regex(pat) {
                    if let Some(caps) = re.captures(bam1) {
                        if let Some(m) = caps.get(1) {
                            base1 = m.as_str().to_string();
                        }
                    }
                }
            }
            base1
        } else if let Some(exp) = explicit_sample {
            exp.to_string()
        } else if emit_intermediate_or_splice {
            extract_sample_name_clean(bam_raw, clean_pat)
        } else {
            format!(
                "{}|{}",
                extract_sample_name_clean(bam1, clean_pat),
                extract_sample_name_clean(bam2, clean_pat)
            )
        }
    } else if let Some(exp) = explicit_sample {
        exp.to_string()
    } else {
        extract_sample_name_clean(bam_raw, clean_pat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_vardict() {
        let c = Config::default();
        assert!((c.freq - 0.01).abs() < f64::EPSILON);
        assert_eq!(c.min_reads, 2);
        assert!((c.base_qual - 22.5).abs() < f64::EPSILON);
        assert_eq!(c.realign, 1);
        assert_eq!(c.flag_filter, 0x504);
        assert_eq!(c.samfilter_raw, "0x504");
        assert!((c.mfreq - 0.25).abs() < f64::EPSILON);
        assert!((c.nmfreq - 0.1).abs() < f64::EPSILON);
        assert_eq!(c.col_chr, 2);
        assert_eq!(c.col_start, 6);
        assert_eq!(c.col_end, 7);
        assert_eq!(c.col_seg_start, 9);
        assert_eq!(c.col_seg_end, 10);
        assert_eq!(c.col_gene, Some(12));
        assert_eq!(c.column_for_chromosome, None);
        assert_eq!(c.trim_bases_after, 0);
        assert_eq!(c.sample_regexp, None);
        assert!(!c.index_csi);
        assert_eq!(c.compat, COMPAT_VARDICTJAVA_1_8_3);
        assert!(!c.vcf_all_alleles);
        assert!(c.af_threshold_raw.is_none());
        assert!(!c.amplicon_single_primer);
        assert_eq!(c.vcf_min_splitreads, 1);
        assert!(!c.vcf_min_splitreads_explicit);

        let core_desc = env!("CARGO_PKG_DESCRIPTION");
        assert!(
            core_desc.is_ascii(),
            "vardictrust-core description must be pure ASCII: {core_desc}"
        );
        assert!(
            core_desc.chars().count() <= 80,
            "vardictrust-core description must be <= 80 chars (got {}): {core_desc}",
            core_desc.chars().count()
        );
        assert!(
            !core_desc.ends_with('.'),
            "vardictrust-core description must not end with '.': {core_desc}"
        );
        assert!(
            core_desc.contains("VarDictRust"),
            "vardictrust-core description must contain 'VarDictRust': {core_desc}"
        );
        assert!(
            !core_desc.contains("--compat"),
            "vardictrust-core description must not contain '--compat': {core_desc}"
        );
    }

    #[test]
    fn compat_wave2_switches_and_min_splitreads_explicit_contract() {
        let vj = std::hint::black_box(COMPAT_VARDICTJAVA_1_8_3);
        assert!(!vj.cigar_homopolymer_merge);
        assert!(!vj.eqx_cigar_read_length);
        assert!(!vj.preserve_discordant_svs);
        assert!(!vj.mark_dup_sv_inclusive_end);
        assert!(!vj.sv_lazy_mate_pileup);
        assert!(!vj.exact_binary_half_even);

        let off = std::hint::black_box(COMPAT_OFF);
        assert!(off.cigar_homopolymer_merge);
        assert!(off.eqx_cigar_read_length);
        assert!(off.preserve_discordant_svs);
        assert!(off.mark_dup_sv_inclusive_end);
        assert!(off.sv_lazy_mate_pileup);
        assert!(off.exact_binary_half_even);

        let mut cfg = Config::default();
        assert!(!cfg.vcf_min_splitreads_explicit);
        cfg.vcf_min_splitreads_explicit = true;
        assert!(cfg.vcf_min_splitreads_explicit);
    }

    #[test]
    fn compat_target_and_profile_contract() {
        assert_eq!(CompatTarget::default(), CompatTarget::VardictJava183);
        assert_eq!(CompatTarget::VardictJava183.to_string(), "vardictjava");
        assert_eq!(CompatTarget::Off.to_string(), "off");
        assert_eq!(
            CompatTarget::VardictJava183.profile(),
            COMPAT_VARDICTJAVA_1_8_3
        );
        assert_eq!(CompatTarget::Off.profile(), COMPAT_OFF);
        assert_eq!(CompatProfile::default(), COMPAT_VARDICTJAVA_1_8_3);

        let vj = std::hint::black_box(COMPAT_VARDICTJAVA_1_8_3);
        let CompatProfile {
            target: vj_target,
            sanitize_sample_info,
            standard_dup_symbolic_alt,
            shift3_3prime_alleles,
            strict_region_bounds,
            suppress_pileup_sentinels,
            live_duprate,
            preserve_mixed_sex_contigs,
            cigar_homopolymer_merge,
            eqx_cigar_read_length,
            preserve_discordant_svs,
            mark_dup_sv_inclusive_end,
            sv_lazy_mate_pileup,
            exact_binary_half_even,
            reconcile_vd_dp,
            coordinate_sorted_vars,
            exact_35_overlap,
            exact_cigar_dlen_match,
            somatic_continue_after_bad_var,
            adaptor_all_consensus_dirs,
            canonical_sv_vcf_anchoring,
            paired_somatic_amplicon,
            cigar_strip_outer_hard_clips,
            nonnegative_subtracted_means,
            reconcile_strand_and_hq_sums,
            reconcile_rd_vd_dp,
            sanitize_iupac_alt_alleles,
            sv_respect_min_len,
            pileup_preserve_shifted_ref,
            recompute_vartype_after_adj_complex,
            strict_vcf_ad_number_r,
            omit_nonsv_vcf_end,
            canonical_vcf_gt,
        } = vj;
        assert_eq!(vj_target, CompatTarget::VardictJava183);
        let vj_switches: [bool; 32] = [
            sanitize_sample_info,
            standard_dup_symbolic_alt,
            shift3_3prime_alleles,
            strict_region_bounds,
            suppress_pileup_sentinels,
            live_duprate,
            preserve_mixed_sex_contigs,
            cigar_homopolymer_merge,
            eqx_cigar_read_length,
            preserve_discordant_svs,
            mark_dup_sv_inclusive_end,
            sv_lazy_mate_pileup,
            exact_binary_half_even,
            reconcile_vd_dp,
            coordinate_sorted_vars,
            exact_35_overlap,
            exact_cigar_dlen_match,
            somatic_continue_after_bad_var,
            adaptor_all_consensus_dirs,
            canonical_sv_vcf_anchoring,
            paired_somatic_amplicon,
            cigar_strip_outer_hard_clips,
            nonnegative_subtracted_means,
            reconcile_strand_and_hq_sums,
            reconcile_rd_vd_dp,
            sanitize_iupac_alt_alleles,
            sv_respect_min_len,
            pileup_preserve_shifted_ref,
            recompute_vartype_after_adj_complex,
            strict_vcf_ad_number_r,
            omit_nonsv_vcf_end,
            canonical_vcf_gt,
        ];
        assert!(vj_switches.iter().all(|&flag| !flag));

        let off = std::hint::black_box(COMPAT_OFF);
        let CompatProfile {
            target: off_target,
            sanitize_sample_info,
            standard_dup_symbolic_alt,
            shift3_3prime_alleles,
            strict_region_bounds,
            suppress_pileup_sentinels,
            live_duprate,
            preserve_mixed_sex_contigs,
            cigar_homopolymer_merge,
            eqx_cigar_read_length,
            preserve_discordant_svs,
            mark_dup_sv_inclusive_end,
            sv_lazy_mate_pileup,
            exact_binary_half_even,
            reconcile_vd_dp,
            coordinate_sorted_vars,
            exact_35_overlap,
            exact_cigar_dlen_match,
            somatic_continue_after_bad_var,
            adaptor_all_consensus_dirs,
            canonical_sv_vcf_anchoring,
            paired_somatic_amplicon,
            cigar_strip_outer_hard_clips,
            nonnegative_subtracted_means,
            reconcile_strand_and_hq_sums,
            reconcile_rd_vd_dp,
            sanitize_iupac_alt_alleles,
            sv_respect_min_len,
            pileup_preserve_shifted_ref,
            recompute_vartype_after_adj_complex,
            strict_vcf_ad_number_r,
            omit_nonsv_vcf_end,
            canonical_vcf_gt,
        } = off;
        assert_eq!(off_target, CompatTarget::Off);
        let off_switches: [bool; 32] = [
            sanitize_sample_info,
            standard_dup_symbolic_alt,
            shift3_3prime_alleles,
            strict_region_bounds,
            suppress_pileup_sentinels,
            live_duprate,
            preserve_mixed_sex_contigs,
            cigar_homopolymer_merge,
            eqx_cigar_read_length,
            preserve_discordant_svs,
            mark_dup_sv_inclusive_end,
            sv_lazy_mate_pileup,
            exact_binary_half_even,
            reconcile_vd_dp,
            coordinate_sorted_vars,
            exact_35_overlap,
            exact_cigar_dlen_match,
            somatic_continue_after_bad_var,
            adaptor_all_consensus_dirs,
            canonical_sv_vcf_anchoring,
            paired_somatic_amplicon,
            cigar_strip_outer_hard_clips,
            nonnegative_subtracted_means,
            reconcile_strand_and_hq_sums,
            reconcile_rd_vd_dp,
            sanitize_iupac_alt_alleles,
            sv_respect_min_len,
            pileup_preserve_shifted_ref,
            recompute_vartype_after_adj_complex,
            strict_vcf_ad_number_r,
            omit_nonsv_vcf_end,
            canonical_vcf_gt,
        ];
        assert!(off_switches.iter().all(|&flag| flag));

        for valid_vj in [
            "vardictjava",
            "VARDICTJAVA",
            "vardictjava-1.8.3",
            "VarDictJava-1.8.3",
            "1.8.3",
            " vardictjava ",
        ] {
            assert_eq!(
                parse_compat_target(valid_vj).unwrap(),
                CompatTarget::VardictJava183,
                "expected VardictJava183 for {valid_vj:?}"
            );
        }
        for valid_off in ["off", "OFF", " off ", "\toff\n"] {
            assert_eq!(
                parse_compat_target(valid_off).unwrap(),
                CompatTarget::Off,
                "expected Off for {valid_off:?}"
            );
        }
        for bad in ["modern", "foo", "", "-G"] {
            let err = parse_compat_target(bad).expect_err("expected Unsupported for bad compat");
            assert!(
                matches!(err, crate::types::VardictError::Unsupported(_)),
                "expected VardictError::Unsupported for {bad:?}, got {err:?}"
            );
            let expected_msg = format!(
                "--compat {bad} not supported; expected 'vardictjava' (or 'vardictjava-1.8.3' / '1.8.3') or 'off'"
            );
            assert_eq!(err.to_string(), format!("unsupported: {expected_msg}"));
        }

        // Cross-module Wave-1/Wave-2 contract checks
        assert_eq!(
            crate::output::vcf::reorder_chromosomes_with_compat(&["chrX", "X", "chr1"], false),
            crate::output::vcf::perl_reorder_chromosomes(&["chrX", "X", "chr1"])
        );
        assert_eq!(
            crate::output::vcf::reorder_chromosomes_with_compat(&["chrX", "X", "chr1"], true),
            vec!["chr1", "X", "chrX"]
        );
        let reference = crate::reference::Reference::new();
        let region = crate::regions::Region {
            chr: "chr1".to_string(),
            start: crate::types::Pos1(1),
            end: crate::types::Pos1(10),
            gene: "G1".to_string(),
            amplicon: None,
        };
        let cfg = Config::default();
        let mut parser =
            crate::cigar::parse::CigarParser::new(&reference, &region, &cfg, false, 0, 0, 0);
        parser.set_live_read_counts(100, 5);
    }

    #[test]
    fn java_integer_decode_table() {
        assert_eq!(java_integer_decode("0x504").unwrap(), 0x504);
        assert_eq!(java_integer_decode("0X504").unwrap(), 0x504);
        assert_eq!(java_integer_decode("#504").unwrap(), 0x504);
        assert_eq!(java_integer_decode("1284").unwrap(), 1284);
        assert_eq!(java_integer_decode("+1284").unwrap(), 1284);
        assert_eq!(java_integer_decode("02404").unwrap(), 0o2404);
        assert_eq!(java_integer_decode("0").unwrap(), 0);
        assert_eq!(java_integer_decode("00").unwrap(), 0);
        assert_eq!(java_integer_decode("-0x10").unwrap(), -16);
        assert_eq!(java_integer_decode("-010").unwrap(), -8);
        assert!(java_integer_decode("0x-1").is_err());
        assert!(java_integer_decode("0x+1").is_err());
        assert!(java_integer_decode("08").is_err());
        assert!(java_integer_decode("0xFFFFFFFF").is_err());
        assert!(java_integer_decode(" 10").is_err());
        assert!(java_integer_decode("").is_err());
    }

    #[test]
    fn extract_sample_name_and_pileup_overrides() {
        assert_eq!(
            extract_sample_name("/path/to/MySample.sorted.bam", None),
            "MySample"
        );
        assert_eq!(
            extract_sample_name("sample_ABC_123.bam", Some(r"sample_([A-Z]+)_\d+")),
            "ABC"
        );
        assert_eq!(extract_sample_name("", None), "null");

        let mut cfg = Config::default();
        cfg.apply_pileup_overrides();
        assert!((cfg.freq - 0.01).abs() < f64::EPSILON);
        assert_eq!(cfg.min_reads, 2);

        cfg.pileup = true;
        cfg.apply_pileup_overrides();
        assert!((cfg.freq - (-1.0)).abs() < f64::EPSILON);
        assert_eq!(cfg.min_reads, 0);
    }

    #[test]
    fn extract_sample_name_java_exact() {
        assert_eq!(
            extract_sample_name("/dir/foo_bar.sorted.dup.bam", None),
            "bar"
        );
        assert_eq!(extract_sample_name("/dir/foo.sorted.dup.bam", None), "foo");
        assert_eq!(
            extract_sample_name("/dir/fooXsorted.bam", None),
            "foo",
            "SAMPLE_PATTERN has unescaped '.' before 'sorted'"
        );
        assert_eq!(
            extract_sample_name("aligned_bams/NA12878.bam", None),
            "aligned",
            "SAMPLE_PATTERN2 matches 'aligned_bam' in directory name when applied to raw -b path"
        );
        assert_eq!(extract_sample_name("/dir/foo_bar.bam", None), "foo_bar");
        assert_eq!(extract_sample_name("/dir/foo.cram", None), "null");
        assert_eq!(
            extract_sample_name("t.bam|n.bam", None),
            "t.bam|n",
            "SAMPLE_PATTERN2 is greedy across `|` when applied to full bamRaw"
        );
        assert_eq!(
            extract_sample_name("/dir/sample_42.bam", Some("/sample_(\\d+)/")),
            "42"
        );
        // -n "/" and -n "//" strip outer slashes first to "" and fall through to SAMPLE_RE1
        assert_eq!(
            extract_sample_name("/dir/foo.sorted.dup.bam", Some("/")),
            "foo"
        );
        assert_eq!(
            extract_sample_name("/dir/foo.sorted.dup.bam", Some("//")),
            "foo"
        );

        // validate_sample_regexp checks
        assert_eq!(validate_sample_regexp(""), Ok(None));
        assert_eq!(validate_sample_regexp("/"), Ok(None));
        assert_eq!(validate_sample_regexp("//"), Ok(None));
        assert_eq!(
            validate_sample_regexp("/sample_(\\d+)/"),
            Ok(Some("sample_(\\d+)".to_string()))
        );
        assert!(validate_sample_regexp("([unclosed").is_err());
        assert!(validate_sample_regexp(r"(?<=/)([^/]+)\.bam$").is_err());

        // resolve_sample_names matrix (paired -R .var/-i, paired BED .var/-i, paired VCF, single)
        assert_eq!(
            resolve_sample_names("t.sorted.bam|n.sorted.bam", None, None, true, true),
            "t"
        );
        assert_eq!(
            resolve_sample_names("t.bam|n.bam", None, None, true, true),
            "t.bam|n"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam:t_2.bam|n_1.bam",
                None,
                Some(r":t_(\d+)\.bam"),
                true,
                true
            ),
            "2"
        );
        assert_eq!(
            resolve_sample_names("t_1.bam|n_2.bam", Some("EXP_T|EXP_N"), None, true, true),
            "EXP_T"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam|n_2.bam",
                Some("EXP_T|EXP_N"),
                Some(r"t_(\d+)"),
                true,
                true
            ),
            "1"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam|n_2.bam",
                Some("EXP_T|EXP_N"),
                Some(r"t_(\d+)"),
                false,
                true
            ),
            "EXP_T|EXP_N"
        );
        assert_eq!(
            resolve_sample_names("t.bam|n.bam", None, None, false, true),
            "t.bam|n"
        );
        assert_eq!(
            resolve_sample_names(
                "t_1.bam|n_2.bam",
                Some("EXP_T|EXP_N"),
                Some(r"t_(\d+)"),
                true,
                false
            ),
            "EXP_T|EXP_N"
        );
        assert_eq!(
            resolve_sample_names(
                "dir/T1.sorted.bam|dir/N1.sorted.bam",
                None,
                None,
                false,
                false
            ),
            "T1|N1"
        );
        assert_eq!(
            resolve_sample_names(
                "dir/s_T9.bam|dir/s_N8.bam",
                None,
                Some("s_([A-Z0-9]+)\\.bam"),
                true,
                false
            ),
            "T9|N8"
        );
        assert_eq!(
            resolve_sample_names(
                "dir/s_T9.bam",
                Some("EXPLICIT"),
                Some("s_([A-Z0-9]+)"),
                false,
                false
            ),
            "EXPLICIT"
        );
        assert_eq!(
            resolve_sample_names("a_1.bam:b_2.bam", None, Some(r":b_(\d+)\.bam"), false, true),
            "2"
        );
        // Single outer-slash strip (`CmdParser.java:86-93`): escaped trailing `\/` and `//...//` anchors
        assert_eq!(
            extract_sample_name("/mnt/cohort/SAMP_ESC/aln.bam", Some(r"/cohort\/([^/]+)\//")),
            "SAMP_ESC"
        );
        assert_eq!(
            resolve_sample_names(
                "/mnt/cohort/T_ESC/t.bam|/mnt/cohort/N_ESC/n.bam",
                None,
                Some(r"/cohort\/([^/]+)\//"),
                false,
                false,
            ),
            "T_ESC|N_ESC"
        );
        assert_eq!(
            resolve_sample_names(
                "/mnt/cohort/T_ESC/t.bam|/mnt/cohort/N_ESC/n.bam",
                None,
                Some(r"/cohort\/([^/]+)\//"),
                true,
                true,
            ),
            "T_ESC"
        );
        assert_eq!(
            extract_sample_name("/cohort/ANCHORED/aln.bam", Some(r"//cohort/([^/]+)//")),
            "ANCHORED"
        );
        assert_eq!(
            extract_sample_name("not_cohort/FAKE/aln_1.bam", Some(r"//cohort/([^/]+)//")),
            "aln_1",
            "double-slash pattern strips once to `/cohort/([^/]+)/` and must not match `not_cohort/`"
        );
    }

    #[test]
    fn vcf_threshold_constructors_and_defaults() {
        let c = Config::default();
        assert!(c.vcf_af.is_none());
        assert!(c.vcf_min_depth.is_none());
        assert!(c.vcf_min_vd.is_none());
        assert!(c.vcf_qmean.is_none());
        assert!(c.vcf_mapq.is_none());
        assert!(c.vcf_pmean.is_none());
        assert!(c.vcf_sn.is_none());
        assert!(c.vcf_max_nm.is_none());
        assert!(c.vcf_gt_freq.is_none());
        assert!(c.vcf_msi.is_none());
        assert!(c.vcf_cluster_dist.is_none());

        let t_new = VcfThreshold::new(20.0_f64);
        assert_eq!(t_new.value, 20.0);
        assert_eq!(t_new.raw, None);

        let t_from: VcfThreshold<f64> = 20.0.into();
        assert_eq!(t_from, t_new);

        let t_raw = VcfThreshold::with_raw(20.0, "20.0".to_string());
        assert_eq!(t_raw.value, 20.0);
        assert_eq!(t_raw.raw.as_deref(), Some("20.0"));
    }

    #[test]
    fn validate_vcf_nonneg_float_token_and_parsers() {
        let (v0, r0) = validate_vcf_nonneg_float_token("+0.00").unwrap();
        assert_eq!(v0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(r0, "0.00");

        let (vn0, rn0) = validate_vcf_nonneg_float_token("-0").unwrap();
        assert_eq!(vn0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rn0, "0");

        let (vn00, rn00) = validate_vcf_nonneg_float_token("-0.0").unwrap();
        assert_eq!(vn00.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rn00, "0.0");

        let (v_sci, r_sci) = validate_vcf_nonneg_float_token("1e-5").unwrap();
        assert_eq!(v_sci.to_bits(), 1e-5_f64.to_bits());
        assert_eq!(r_sci, "1e-5");

        let exact_64 = format!("0.{}", "0".repeat(62));
        assert_eq!(exact_64.len(), 64);
        assert!(validate_vcf_nonneg_float_token(&exact_64).is_ok());

        let over_64 = format!("0.{}", "0".repeat(63));
        assert_eq!(over_64.len(), 65);
        assert!(validate_vcf_nonneg_float_token(&over_64).is_err());

        for bad in [
            "", " 0.05", "0.05 ", "0.05\n", "0.05\t", "0.05;bad", "-0.01", "-1e-400", "-3", "NaN",
            "inf", "Infinity", "1e500", "0.01f", "-", "+", ".",
        ] {
            assert!(
                validate_vcf_nonneg_float_token(bad).is_err(),
                "expected Err for {bad:?}"
            );
        }

        let thr = parse_vcf_threshold_nonneg_f64("+0.00", "--vcf-af").unwrap();
        assert_eq!(thr.value.to_bits(), 0.0_f64.to_bits());
        assert_eq!(thr.raw.as_deref(), Some("0.00"));

        assert_eq!(
            parse_vcf_unit_f64("-0.0", "--vcf-gt-freq")
                .unwrap()
                .to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(parse_vcf_unit_f64("1.0", "--vcf-gt-freq").unwrap(), 1.0);
        assert!(parse_vcf_unit_f64("1.0001", "--vcf-gt-freq").is_err());
        assert!(parse_vcf_unit_f64("-0.1", "--vcf-gt-freq").is_err());

        assert_eq!(parse_vcf_nonneg_i32("0", "--min-splitreads").unwrap(), 0);
        assert_eq!(parse_vcf_nonneg_i32("15", "--min-splitreads").unwrap(), 15);
        for bad_int in ["", " 1", "-1", "-3", "2.5", "0x2"] {
            assert!(
                parse_vcf_nonneg_i32(bad_int, "--min-splitreads").is_err(),
                "expected Err for {bad_int:?}"
            );
            assert!(
                validate_vcf_nonneg_i32_token(bad_int).is_err(),
                "expected validate_vcf_nonneg_i32_token Err for {bad_int:?}"
            );
        }
        assert_eq!(validate_vcf_nonneg_i32_token("03").unwrap(), (3, "03"));
        assert_eq!(validate_vcf_nonneg_i32_token("+02").unwrap(), (2, "+02"));
        assert_eq!(validate_vcf_nonneg_i32_token("010").unwrap(), (10, "010"));
        assert_eq!(validate_vcf_nonneg_i32_token("00").unwrap(), (0, "00"));
        assert_eq!(validate_vcf_nonneg_i32_token("+0").unwrap(), (0, "+0"));
        assert_eq!(validate_vcf_nonneg_i32_token("-0").unwrap(), (0, "0"));
        assert!(validate_vcf_nonneg_i32_token(&"0".repeat(64)).is_ok());
        assert!(validate_vcf_nonneg_i32_token(&"0".repeat(65)).is_err());
        assert_eq!(
            parse_vcf_threshold_nonneg_i32("03", "--vcf-min-depth").unwrap(),
            VcfThreshold::with_raw(3, "03".to_string())
        );
        assert_eq!(
            parse_vcf_threshold_nonneg_i32("+02", "--vcf-min-vd").unwrap(),
            VcfThreshold::with_raw(2, "+02".to_string())
        );
        assert_eq!(
            parse_vcf_threshold_nonneg_i32("010", "--vcf-cluster-dist").unwrap(),
            VcfThreshold::with_raw(10, "010".to_string())
        );
        assert_eq!(
            parse_vcf_threshold_nonneg_i32("-0", "--vcf-cluster-dist").unwrap(),
            VcfThreshold::with_raw(0, "0".to_string())
        );
    }

    #[test]
    fn config_validate_vcf_thresholds_checks_all_fields() {
        let mut cfg = Config {
            af_threshold: -0.0,
            vcf_pval: -0.0,
            vcf_pval_raw: Some("-0.0".to_string()),
            vcf_af: Some(VcfThreshold::with_raw(-0.0, "+0.00".to_string())),
            vcf_min_depth: Some(VcfThreshold::with_raw(3, "03".to_string())),
            vcf_min_vd: Some(VcfThreshold::with_raw(2, "+02".to_string())),
            vcf_qmean: Some(VcfThreshold::with_raw(20.0, "20.0".to_string())),
            vcf_mapq: Some(VcfThreshold::with_raw(0.0, "-0".to_string())),
            vcf_pmean: Some(VcfThreshold::with_raw(7.5, "7.5".to_string())),
            vcf_sn: Some(VcfThreshold::new(2.0)),
            vcf_max_nm: Some(VcfThreshold::with_raw(4.0, "4.0".to_string())),
            vcf_gt_freq: Some(-0.0),
            vcf_msi: Some(VcfThreshold::with_raw(10.5, "10.5".to_string())),
            vcf_cluster_dist: Some(VcfThreshold::with_raw(0, "-0".to_string())),
            ..Config::default()
        };
        cfg.validate_vcf_thresholds().unwrap();
        assert_eq!(cfg.af_threshold.to_bits(), 0.0_f64.to_bits());
        assert_eq!(cfg.vcf_pval.to_bits(), 0.0_f64.to_bits());
        assert_eq!(cfg.vcf_pval_raw.as_deref(), Some("0.0"));
        assert_eq!(
            cfg.vcf_af.as_ref().unwrap().value.to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(cfg.vcf_af.as_ref().unwrap().raw.as_deref(), Some("0.00"));
        assert_eq!(
            cfg.vcf_min_depth.as_ref().unwrap().raw.as_deref(),
            Some("03")
        );
        assert_eq!(cfg.vcf_min_vd.as_ref().unwrap().raw.as_deref(), Some("+02"));
        assert_eq!(
            cfg.vcf_cluster_dist.as_ref().unwrap().raw.as_deref(),
            Some("0")
        );
        assert_eq!(cfg.vcf_mapq.as_ref().unwrap().raw.as_deref(), Some("0"));
        assert_eq!(cfg.vcf_gt_freq.unwrap().to_bits(), 0.0_f64.to_bits());

        // Reject negative integer overrides, >1.0 gt_freq, mismatched/hostile raw tokens
        let mut bad_depth = Config {
            vcf_min_depth: Some(VcfThreshold::new(-1)),
            ..Config::default()
        };
        assert!(bad_depth.validate_vcf_thresholds().is_err());

        let mut bad_vd = Config {
            vcf_min_vd: Some(VcfThreshold::new(-2)),
            ..Config::default()
        };
        assert!(bad_vd.validate_vcf_thresholds().is_err());

        let mut bad_cluster = Config {
            vcf_cluster_dist: Some(VcfThreshold::new(-5)),
            ..Config::default()
        };
        assert!(bad_cluster.validate_vcf_thresholds().is_err());

        let mut bad_int_raw_mismatch = Config {
            vcf_cluster_dist: Some(VcfThreshold::with_raw(10, "015".to_string())),
            ..Config::default()
        };
        assert!(bad_int_raw_mismatch.validate_vcf_thresholds().is_err());

        let mut bad_gt = Config {
            vcf_gt_freq: Some(1.1),
            ..Config::default()
        };
        assert!(bad_gt.validate_vcf_thresholds().is_err());

        let mut bad_raw_mismatch = Config {
            vcf_qmean: Some(VcfThreshold::with_raw(20.0, "25.0".to_string())),
            ..Config::default()
        };
        assert!(bad_raw_mismatch.validate_vcf_thresholds().is_err());

        let mut bad_raw_inject = Config {
            vcf_qmean: Some(VcfThreshold::with_raw(20.0, "20.0\n##evil".to_string())),
            ..Config::default()
        };
        assert!(bad_raw_inject.validate_vcf_thresholds().is_err());

        let mut bad_raw_over64 = Config {
            vcf_msi: Some(VcfThreshold::with_raw(0.0, format!("0.{}", "0".repeat(63)))),
            ..Config::default()
        };
        assert!(bad_raw_over64.validate_vcf_thresholds().is_err());
    }

    #[test]
    fn batch12_parse_java_number_usize_enforces_java_int_bounds() {
        assert_eq!(parse_java_number_usize("0").unwrap(), 0);
        assert_eq!(
            parse_java_number_usize("2147483647").unwrap(),
            i32::MAX as usize
        );
        for bad in [
            "-1",
            "2147483648",
            "18446744073709551615",
            "\x00all_cores",
            "",
            " 1",
        ] {
            assert!(
                parse_java_number_usize(bad).is_err(),
                "expected Err from parse_java_number_usize({bad:?})"
            );
        }
    }

    #[test]
    fn wave4_schema_compat_21_java_float_and_amplicon_concordance_contract() {
        use std::collections::HashSet;

        let amp_a = AmpliconConcordance {
            gdamp: 3,
            tlamp: 4,
            ncamp: 1,
            ampflag: true,
        };
        let amp_b = amp_a;
        assert_eq!(amp_a, amp_b);
        let mut set = HashSet::new();
        assert!(set.insert(amp_a));
        assert!(!set.insert(amp_b));

        for (raw, expected) in [
            ("1e-5", 1e-5),
            ("1E-5", 1e-5),
            ("1e1", 10.0),
            ("1.5e2", 150.0),
            ("0.05d", 0.05),
            ("0.05f", 0.05),
            (" 0.05D ", 0.05),
            (" 1e-5f ", 1e-5),
            (".01", 0.01),
            ("0e0", 0.0),
        ] {
            let got = parse_java_number_f64(raw)
                .unwrap_or_else(|e| panic!("parse_java_number_f64({raw:?}) failed: {e}"));
            assert!(
                (got - expected).abs() <= 1e-15,
                "raw={raw:?}: got {got}, expected {expected}"
            );
        }
        assert!(parse_java_number_i32("1e1").is_err());
        for bad in [
            "1e", "e", "1e500", "-1e500", "NaN", "inf", "Infinity", " 10", "10 ", "",
        ] {
            assert!(
                parse_java_number_f64(bad).is_err(),
                "expected Err from parse_java_number_f64({bad:?})"
            );
        }

        let var_new = crate::variations::Variation::new();
        assert!(!var_new.pre_anchored_dup);
        assert!(!crate::variations::Variation::default().pre_anchored_dup);

        let som = crate::postprocess::somatic::SomaticVariant {
            variant1: crate::variations::Variant::new(),
            variant2: None,
            status: "StrongSomatic".to_string(),
            sv1: String::new(),
            sv2: String::new(),
            end_is_variant2: false,
            tumor_amplicon: Some(amp_a),
            normal_amplicon: None,
        };
        assert_eq!(som.tumor_amplicon, Some(amp_a));
        assert!(som.normal_amplicon.is_none());

        let mut cfg_af = Config {
            af_threshold: -0.0,
            af_threshold_raw: Some("-0.0".to_string()),
            ..Config::default()
        };
        cfg_af.validate_vcf_thresholds().unwrap();
        assert_eq!(cfg_af.af_threshold.to_bits(), 0.0_f64.to_bits());
        assert_eq!(cfg_af.af_threshold_raw.as_deref(), Some("0.0"));

        let mut bad_af_mismatch = Config {
            af_threshold: 0.05,
            af_threshold_raw: Some("0.01".to_string()),
            ..Config::default()
        };
        assert!(bad_af_mismatch.validate_vcf_thresholds().is_err());

        let mut bad_af_inject = Config {
            af_threshold: 0.05,
            af_threshold_raw: Some("0.05;evil".to_string()),
            ..Config::default()
        };
        assert!(bad_af_inject.validate_vcf_thresholds().is_err());
    }
}
