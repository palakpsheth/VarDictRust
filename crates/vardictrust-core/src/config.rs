//! Configuration — mirrors VarDictJava's `Configuration` (§4.2 defaults).

/// Central configuration struct.
///
/// All VarDictJava defaults from §4.2 are preserved here. Constructed from
/// the CLI args via `From<Cli>`.
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

    // -- Additional flags added in T027-T032 --
    /// Extension for BED region (default: 0).
    pub extension: i32,
    /// Number of bases to extend for reads (default: 2).
    pub extension_nuc: i32,
    /// Include 3' UTR.
    pub three_prime: bool,
    /// Minimum # of reads to call with (default: 2).
    pub min_reads_bias: u32,
    /// qratio threshold (default: 1.5).
    pub qratio: f64,
    /// MapQ for output (default: 0).
    pub out_mapq: u32,
    /// Indel size (default: 50).
    pub indel_size: u32,
    /// Zero-based BED coordinates.
    pub zero_based: bool,
    /// SAM flag filter (default: 0x504).
    pub flag_filter: u32,
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
    /// unique mode alignment enabled (-u)
    pub unique_mode_alignment_enabled: bool,
    /// unique mode second in pair enabled (-UN)
    pub unique_mode_second_in_pair_enabled: bool,
    /// Omit END tag in VCF output (matches var2vcf_valid.pl -E)
    pub omit_end_tag: bool,
    /// Include Ns in depth calculation (-K)
    pub include_n: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            exception_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
            reference: String::new(),
            bam: String::new(),
            freq: 0.01,
            af_threshold: 0.01,
            sample: "SAMPLE".to_string(),
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
            extension: 0,
            extension_nuc: 2,
            three_prime: false,
            min_reads_bias: 2,
            qratio: 1.5,
            out_mapq: 0,
            indel_size: 50,
            zero_based: false,
            flag_filter: 0x504,
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
            include_n: false,
        }
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
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
        assert!((c.mfreq - 0.25).abs() < f64::EPSILON);
        assert!((c.nmfreq - 0.1).abs() < f64::EPSILON);
        assert_eq!(c.col_chr, 2);
        assert_eq!(c.col_start, 6);
        assert_eq!(c.col_end, 7);
        assert_eq!(c.col_seg_start, 9);
        assert_eq!(c.col_seg_end, 10);
        assert_eq!(c.col_gene, Some(12));
        assert_eq!(c.column_for_chromosome, None);
    }
}
