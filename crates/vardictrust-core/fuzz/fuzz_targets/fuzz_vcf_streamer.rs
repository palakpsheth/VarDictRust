#![no_main]

use std::collections::HashMap;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::{Config, VcfThreshold};
use vardictrust_core::output::table::{
    amplicon_38_into, amplicon_40_into, simple_36_into, simple_38_into, somatic_55_into,
    somatic_61_into,
};
use vardictrust_core::output::vcf::{
    build_single_vcf_event, effective_vcf_cluster_dist, format_paired_vcf_record,
    format_vcf_record_into, PairedRecordStreamer, SingleRecordStreamer, SingleVcfEvent,
    SingleVcfRecord,
};
use vardictrust_core::postprocess::amplicon::AmpliconVariant;
use vardictrust_core::postprocess::somatic::SomaticVariant;
use vardictrust_core::regions::Region;
use vardictrust_core::types::Pos1;
use vardictrust_core::variations::Variant;

/// Maximum accepted libFuzzer input payload (`8 KiB`) to bound per-iteration runtime.
const MAX_FUZZ_INPUT_BYTES: usize = 8192;
/// Maximum synthesized variant rows per input (`32`).
const MAX_FUZZ_ROWS: usize = 32;
/// Maximum synthesized `REF`/`ALT` allele byte length (`256`).
const MAX_ALLELE_BYTES: usize = 256;
/// Upper bound (`2048`) on normal per-row Fisher counts (`support_len <= 4097`).
const MAX_NORMAL_FISHER_COUNT: u16 = 2048;
/// Floor (`> 16_777_216 = MAX_DENSE_FISHER_SUPPORT`) forcing `FisherExact` unimodal mode without `128 MiB` `logdc` allocation.
const UNIMODAL_SUPPORT_FLOOR: i32 = 16_777_220;

struct ByteCursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> ByteCursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    fn u8(&mut self) -> u8 {
        if self.data.is_empty() {
            return 0;
        }
        let idx = self.pos % self.data.len();
        self.pos = self.pos.wrapping_add(1);
        self.data.get(idx).copied().unwrap_or(0)
    }

    fn bool(&mut self) -> bool {
        (self.u8() & 1) != 0
    }

    fn u16(&mut self) -> u16 {
        u16::from_le_bytes([self.u8(), self.u8()])
    }

    fn i32(&mut self) -> i32 {
        i32::from_le_bytes([self.u8(), self.u8(), self.u8(), self.u8()])
    }

    fn u64(&mut self) -> u64 {
        u64::from_le_bytes([
            self.u8(),
            self.u8(),
            self.u8(),
            self.u8(),
            self.u8(),
            self.u8(),
            self.u8(),
            self.u8(),
        ])
    }

    fn metric_f64(&mut self, max_norm: f64) -> f64 {
        match self.u8() & 0x0f {
            0 => 0.0,
            1 => -0.0,
            2 => f64::from_bits(0x3f0a_36e2_eb1c_432d), // 5e-5 HALF_EVEN tie
            3 => f64::from_bits(0x3ed4_f8b5_88e3_68f1), // 5e-6 HALF_EVEN tie
            4 => f64::from_bits(0x3f40_624d_d2f1_a9fc), // 5e-4 HALF_EVEN tie
            5 => f64::NAN,
            6 => f64::INFINITY,
            7 => f64::NEG_INFINITY,
            8 => f64::from_bits(self.u64()),
            _ => (f64::from(self.u16()) / 65535.0) * max_norm,
        }
    }

    /// UTF-8 string up to `max_bytes` bytes (including `\t`/`\r`/`\n`/`\0` when `allow_ctrl` is true,
    /// or sanitized without `\t`/`\r`/`\n` when `allow_ctrl` is false).
    fn utf8_capped(&mut self, max_bytes: usize, allow_ctrl: bool) -> String {
        let target_len = usize::from(self.u16()) % (max_bytes.saturating_add(1)).max(1);
        let mut raw = Vec::with_capacity(target_len);
        for _ in 0..target_len {
            raw.push(self.u8());
        }
        let lossy = String::from_utf8_lossy(&raw);
        let mut out = String::with_capacity(target_len);
        for ch in lossy.chars() {
            let clean = if !allow_ctrl && matches!(ch, '\t' | '\r' | '\n') {
                '_'
            } else {
                ch
            };
            if out.len().saturating_add(clean.len_utf8()) > max_bytes {
                break;
            }
            out.push(clean);
        }
        out
    }

    /// UTF-8 token up to `max_bytes` bytes without `\t`, `\r`, `\n`, or `;` (for non-SAMPLE VCF INFO values).
    fn utf8_no_semi(&mut self, max_bytes: usize) -> String {
        let s = self.utf8_capped(max_bytes, false);
        s.replace(';', "_")
    }

    /// Synthesizes an allele `<= 256` bytes covering empty/0 Perl-skip tokens, SNVs, indels,
    /// symbolic SV alleles (`<DUP>`, `<INV>`, `<DEL>`, `<INS>`, `<dup-N>`), and UTF-8 tokens.
    fn allele_256(&mut self) -> String {
        const BASES: &[u8; 5] = b"ACGTN";
        match self.u8() % 12 {
            0 => String::new(),
            1 => "0".to_string(),
            2 => String::from(BASES[usize::from(self.u8()) % BASES.len()] as char),
            3 => "<DUP>".to_string(),
            4 => "<INV>".to_string(),
            5 => "<DEL>".to_string(),
            6 => "<INS>".to_string(),
            7 => format!("AC<dup-{}>GT", (self.u8() % 20) + 1),
            8..=10 => {
                let len = ((usize::from(self.u8()) % MAX_ALLELE_BYTES) + 1).min(MAX_ALLELE_BYTES);
                let mut s = String::with_capacity(len);
                for _ in 0..len {
                    s.push(BASES[usize::from(self.u8()) % BASES.len()] as char);
                }
                s
            }
            _ => self.utf8_no_semi(MAX_ALLELE_BYTES),
        }
    }

    fn opt_float_threshold(&mut self, default_max: f64) -> Option<VcfThreshold<f64>> {
        match self.u8() % 4 {
            0 => None,
            1 => Some(VcfThreshold::new(self.metric_f64(default_max))),
            2 => {
                let v = f64::from(self.u8() % 50) * 0.5;
                let raw = format!("{v:.2}");
                Some(VcfThreshold::with_raw(v, raw))
            }
            _ => {
                let v = self.metric_f64(default_max);
                let raw = self.utf8_capped(72, true);
                Some(VcfThreshold::with_raw(v, raw))
            }
        }
    }

    fn opt_int_threshold(&mut self, max_v: i32) -> Option<VcfThreshold<i32>> {
        match self.u8() % 4 {
            0 => None,
            1 => Some(VcfThreshold::new(self.i32() % max_v.max(1))),
            2 => {
                let v = i32::from(self.u8() % 32);
                Some(VcfThreshold::with_raw(v, format!("+0{v}")))
            }
            _ => {
                let v = self.i32();
                let raw = self.utf8_capped(72, true);
                Some(VcfThreshold::with_raw(v, raw))
            }
        }
    }
}

#[inline]
fn assert_tsv_cols(kind: &str, line: &str, expected_cols: usize) {
    assert!(
        std::str::from_utf8(line.as_bytes()).is_ok(),
        "{kind} row is not valid UTF-8"
    );
    assert!(
        !line.bytes().any(|b| matches!(b, b'\r' | b'\n')),
        "{kind} row contains embedded CR/LF"
    );
    let actual_cols = line.bytes().filter(|&b| b == b'\t').count() + 1;
    assert_eq!(
        actual_cols, expected_cols,
        "{kind} column count mismatch: got {actual_cols}, expected {expected_cols}"
    );
}

fn assert_vcf_data_record(kind: &str, line: &str, expected_cols: usize, config: &Config) {
    assert_tsv_cols(kind, line, expected_cols);
    if !config.compat.sanitize_sample_info || line.starts_with('#') {
        return;
    }
    let cols: Vec<&str> = line.split('\t').collect();
    let alt = cols.get(4).copied().unwrap_or("");
    let info = cols.get(7).copied().unwrap_or("");
    let entries: Vec<&str> = info.split(';').collect();

    let sample_idx = if expected_cols == 10 { 0usize } else { 1usize };
    let sample_kv = entries.get(sample_idx).copied().unwrap_or("");
    let Some(sample_val) = sample_kv.strip_prefix("SAMPLE=") else {
        panic!("{kind} missing SAMPLE= at INFO index {sample_idx}: {info:?}");
    };
    assert!(
        !sample_val.is_empty() && sample_val != ".",
        "{kind} has empty or '.' INFO/SAMPLE value: {sample_val:?}"
    );
    assert!(
        !sample_val.contains(',')
            && !sample_val.contains('=')
            && !sample_val.contains(';')
            && !sample_val.chars().any(char::is_whitespace),
        "{kind} has unescaped VCF delimiter or whitespace in INFO/SAMPLE: {sample_val:?}"
    );
    assert_eq!(
        entries
            .iter()
            .filter(|kv| kv.starts_with("SAMPLE="))
            .count(),
        1,
        "{kind} must contain exactly one SAMPLE= key in INFO: {info:?}"
    );
    assert!(
        entries
            .get(sample_idx + 1)
            .is_some_and(|kv| kv.starts_with("TYPE=")),
        "{kind} INFO key immediately after SAMPLE= must be TYPE=: {info:?}"
    );

    let expected_keys = if expected_cols == 10 {
        let base_end = 26 + usize::from(!config.omit_end_tag);
        if entries.iter().any(|kv| kv.starts_with("GDAMP=")) {
            base_end + 4
        } else {
            let has_crispr = usize::from(matches!(config.crispr_cutting_site, Some(n) if n != 0));
            let sv_extra = if alt != "." && alt.contains('<') {
                2
            } else {
                0
            };
            base_end + 3 + has_crispr + sv_extra
        }
    } else {
        let flag_extra = usize::from(
            config.vcf_paired_strict
                && matches!(
                    entries.last().copied(),
                    Some("P0.01Likely" | "InDelLikely")
                ),
        );
        13 + flag_extra
    };
    assert_eq!(
        entries.len(),
        expected_keys,
        "{kind} INFO key count mismatch: got {}, expected {expected_keys} in {info:?}",
        entries.len()
    );
}

fn assert_vcf_header(kind: &str, hdr_bytes: &[u8], expected_chrom_cols: usize, compat_off: bool) {
    let hdr_str = std::str::from_utf8(hdr_bytes).expect("VCF header must be valid UTF-8");
    assert!(
        !hdr_str
            .chars()
            .any(|c| matches!(c, '\r' | '\u{0085}' | '\u{2028}' | '\u{2029}')),
        "{kind} contains forbidden CR or Unicode line separator"
    );
    assert!(
        hdr_str.ends_with('\n'),
        "{kind} must end with trailing newline"
    );
    let trimmed = &hdr_str[..hdr_str.len() - 1];
    let mut chrom_count = 0usize;
    let lines: Vec<&str> = trimmed.split('\n').collect();
    assert!(!lines.is_empty(), "{kind} emitted zero lines");
    for (idx, line) in lines.iter().enumerate() {
        if idx + 1 < lines.len() {
            assert!(
                line.starts_with("##"),
                "{kind} meta line {idx} does not start with '##': {line:?}"
            );
        } else {
            assert!(
                line.starts_with("#CHROM\t"),
                "{kind} final header line does not start with '#CHROM\\t': {line:?}"
            );
        }
        if line.starts_with("#CHROM") {
            chrom_count += 1;
            assert_tsv_cols(kind, line, expected_chrom_cols);
        }
    }
    assert_eq!(
        chrom_count, 1,
        "{kind} must contain exactly one #CHROM line, got {chrom_count}"
    );
    if compat_off {
        assert_eq!(
            lines
                .iter()
                .filter(|&&l| l == "##VarDictRustCompat=off")
                .count(),
            1,
            "{kind} must emit ##VarDictRustCompat=off once when compat=off"
        );
        assert_eq!(lines.get(2), Some(&"##VarDictRustCompat=off"));
    } else {
        assert!(
            !lines.iter().any(|l| l.starts_with("##VarDictRustCompat")),
            "{kind} must not emit ##VarDictRustCompat under default compat"
        );
    }
}

fn synth_variant(
    cur: &mut ByteCursor<'_>,
    pos_override: Option<i32>,
    allow_unimodal: bool,
    unimodal_used: &mut bool,
) -> Variant {
    let mut v = Variant::new();
    let count_mod = MAX_NORMAL_FISHER_COUNT + 1;

    if allow_unimodal && !*unimodal_used {
        *unimodal_used = true;
        if cur.bool() {
            // Unimodal via support_len > 16_777_216 (MAX_DENSE_FISHER_SUPPORT).
            let j0 = i32::from(cur.u16() % 512);
            let j1 = i32::from(cur.u16() % 512);
            v.ref_fwd_cov = UNIMODAL_SUPPORT_FLOOR + j0;
            v.rev = UNIMODAL_SUPPORT_FLOOR + j1;
            v.ref_rev_cov = i32::from(cur.u16() % count_mod);
            v.fwd = i32::from(cur.u16() % count_mod);
        } else {
            // Unimodal via m64 + n64 > i32::MAX (!fits_dense_i32) with counts > 16_777_216.
            v.ref_fwd_cov = 1_100_000_000 + i32::from(cur.u16() % 4096);
            v.ref_rev_cov = 1_100_000_000 + i32::from(cur.u16() % 4096);
            v.fwd = i32::from(cur.u16() % 512);
            v.rev = i32::from(cur.u16() % 512);
        }
        v.pos_coverage = v.fwd.saturating_add(v.rev);
        v.total_pos_cov = v
            .pos_coverage
            .saturating_add(v.ref_fwd_cov)
            .saturating_add(v.ref_rev_cov);
    } else {
        // Strictly cap per-row Fisher counts at 0..=2048.
        if (cur.u8() & 0x0f) == 0 {
            v.ref_fwd_cov = 0;
            v.ref_rev_cov = 0;
            v.fwd = 0;
            v.rev = 0;
            v.pos_coverage = 0;
            v.total_pos_cov = 0;
        } else {
            v.ref_fwd_cov = i32::from(cur.u16() % count_mod);
            v.ref_rev_cov = i32::from(cur.u16() % count_mod);
            v.fwd = i32::from(cur.u16() % count_mod);
            v.rev = i32::from(cur.u16() % count_mod);
            v.pos_coverage = i32::from(cur.u16() % count_mod);
            // Keep tref = clamp(total_pos_cov - pos_coverage, 0, ..) in 0..=2048
            // while also exercising VD > DP (total_pos_cov < pos_coverage).
            let delta = i32::from(cur.u16() % (MAX_NORMAL_FISHER_COUNT + 257)) - 256;
            v.total_pos_cov = v.pos_coverage.saturating_add(delta).max(0);
        }
    }

    let start_pos = pos_override.unwrap_or_else(|| match cur.u8() & 7 {
        0 => 0,
        1 => -1,
        2 => i32::MIN,
        3 => i32::MAX,
        _ => i32::from(cur.u16()).saturating_add(1),
    });
    v.start_pos = start_pos;
    let span = i32::from(cur.u8() as i8);
    v.end_pos = start_pos.saturating_add(span);

    v.ref_allele = cur.allele_256();
    v.var_allele = if (cur.u8() & 7) == 0 {
        v.ref_allele.clone()
    } else {
        cur.allele_256()
    };
    debug_assert!(v.ref_allele.len() <= 256 && v.var_allele.len() <= 256);

    v.bias = match cur.u8() % 7 {
        0 => "0".to_string(),
        1 => "0;0".to_string(),
        2 => "2;1".to_string(),
        3 => "2;2".to_string(),
        4 => "1;2".to_string(),
        5 => "0;2".to_string(),
        _ => cur.utf8_capped(12, false),
    };
    v.var_type = match cur.u8() % 9 {
        0 => String::new(),
        1 => "SNV".to_string(),
        2 => "Insertion".to_string(),
        3 => "Deletion".to_string(),
        4 => "Complex".to_string(),
        5 => "MNV".to_string(),
        6 => "DUP".to_string(),
        7 => "INV".to_string(),
        _ => cur.utf8_no_semi(12),
    };
    v.genotype = match cur.u8() % 6 {
        0 => None,
        1 => Some(String::new()),
        2 => Some("0".to_string()),
        3 => Some("1/1".to_string()),
        4 => Some("A/G".to_string()),
        _ => Some(cur.utf8_capped(12, false)),
    };
    v.sv = match cur.u8() % 5 {
        0 => String::new(),
        1 => "0".to_string(),
        2 => format!("{}-{}-{}", cur.u8() % 10, cur.u8() % 20, cur.u8() % 8),
        3 => "1-2-0".to_string(),
        _ => cur.utf8_capped(24, false),
    };
    v.left_seq = match cur.u8() % 4 {
        0 => String::new(),
        1 => "A".repeat(16), // LongAT trigger (>= 14)
        2 => "ACGTACGT".to_string(),
        _ => cur.utf8_no_semi(32),
    };
    v.right_seq = match cur.u8() % 4 {
        0 => String::new(),
        1 => "T".repeat(16),
        2 => "TGCATGCA".to_string(),
        _ => cur.utf8_no_semi(32),
    };

    v.frequency = cur.metric_f64(1.05);
    v.mean_position = cur.metric_f64(50.0);
    v.pstd = cur.bool();
    v.mean_quality = cur.metric_f64(45.0);
    v.qstd = cur.bool();
    v.mean_mapping_quality = cur.metric_f64(60.0);
    v.hq_lq_ratio = cur.metric_f64(20.0);
    v.hi_freq = cur.metric_f64(1.05);
    v.extra_freq = cur.metric_f64(1.05);
    v.shift3 = i32::from(cur.u8() % 32) - 4;
    v.msi = cur.metric_f64(20.0);
    v.msint = i32::from(cur.u8() % 16) - 2;
    v.nm = cur.metric_f64(12.0);
    v.hicnt = i32::from(cur.u16() % 2049);
    v.hicov = i32::from(cur.u16() % 2049);
    v.duprate = cur.metric_f64(1.0);
    v.crispr = cur.i32();
    v.null_record = (cur.u8() & 7) == 0;
    v.debug = cur.utf8_capped(24, false);
    v
}

fuzz_target!(|data: &[u8]| {
    if data.len() < 16 || data.len() > MAX_FUZZ_INPUT_BYTES {
        return;
    }

    let mut cur = ByteCursor::new(data);
    let mode_bits = cur.u8();
    let num_rows = ((usize::from(cur.u8()) % MAX_FUZZ_ROWS) + 1).min(MAX_FUZZ_ROWS);
    let clustered_positions = (mode_bits & 1) == 0;
    let enable_unimodal_once = (mode_bits & 0x3f) == 0x3f;
    let unimodal_target_row = usize::from(cur.u8()) % num_rows;

    let config = Config {
        fisher: (mode_bits & 2) != 0,
        vcf_pass_only: (mode_bits & 4) != 0,
        vcf_paired_strict: (mode_bits & 8) != 0,
        omit_end_tag: (mode_bits & 16) != 0,
        vcf_pstd_filter: (mode_bits & 32) != 0,
        af_threshold_truthy_zero: cur.bool(),
        af_threshold: cur.metric_f64(0.5),
        vcf_min_splitreads: i32::from(cur.u8() % 8) - 1,
        crispr_cutting_site: match cur.u8() % 3 {
            0 => None,
            1 => Some(0),
            _ => Some(cur.i32()),
        },
        sample: match cur.u8() % 8 {
            0 => "SAMPLE".to_string(),
            1 => "TUMOR|NORMAL".to_string(),
            2 => "TUMOR|".to_string(),
            3 => String::new(),
            4 => ".".to_string(),
            5 => "..|.".to_string(),
            6 => "a;FOO=1,bar%20baz\tq|NORM;X=2".to_string(),
            _ => cur.utf8_capped(48, true),
        },
        vcf_pval: cur.metric_f64(0.5),
        vcf_pval_raw: match cur.u8() % 3 {
            0 => None,
            1 => Some("0.050".to_string()),
            _ => Some(cur.utf8_capped(72, true)),
        },
        vcf_gt_freq: if cur.bool() {
            Some(cur.metric_f64(1.2))
        } else {
            None
        },
        vcf_af: cur.opt_float_threshold(1.0),
        vcf_min_depth: cur.opt_int_threshold(100),
        vcf_min_vd: cur.opt_int_threshold(50),
        vcf_qmean: cur.opt_float_threshold(45.0),
        vcf_mapq: cur.opt_float_threshold(60.0),
        vcf_pmean: cur.opt_float_threshold(30.0),
        vcf_sn: cur.opt_float_threshold(10.0),
        vcf_max_nm: cur.opt_float_threshold(12.0),
        vcf_msi: cur.opt_float_threshold(20.0),
        vcf_cluster_dist: cur.opt_int_threshold(50),
        compat: if cur.bool() {
            vardictrust_core::config::COMPAT_OFF
        } else {
            vardictrust_core::config::COMPAT_VARDICTJAVA_1_8_3
        },
        vcf_all_alleles: cur.bool(),
        ..Default::default()
    };

    let _ = config.clone().validate_vcf_thresholds();
    let opt_c = effective_vcf_cluster_dist(&config);

    // Table region may contain control chars to exercise `push_tsv_clean`;
    // VCF region uses sanitized non-tab/newline contig name.
    let table_sample = cur.utf8_capped(32, true);
    let reg_start = (u32::from(cur.u16()) % 10_000).saturating_add(1);
    let reg_end = reg_start.saturating_add(u32::from(cur.u16() % 4096));
    let table_region = Region {
        chr: {
            let s = cur.utf8_capped(24, true);
            if s.is_empty() {
                "chr1".to_string()
            } else {
                s
            }
        },
        start: Pos1(reg_start),
        end: Pos1(reg_end),
        gene: cur.utf8_capped(24, true),
        amplicon: None,
    };
    let vcf_chr = {
        let s = cur.utf8_capped(24, false);
        if s.is_empty() {
            "chr1".to_string()
        } else {
            s
        }
    };
    let vcf_region = Region {
        chr: vcf_chr.clone(),
        start: Pos1(reg_start),
        end: Pos1(reg_end),
        gene: "FUZZ".to_string(),
        amplicon: None,
    };

    let mut single_streamer_cfg = SingleRecordStreamer::from_config(&config);
    let mut single_streamer_raw = SingleRecordStreamer::with_pass_only(opt_c, config.vcf_pass_only);
    let mut paired_streamer_cfg = PairedRecordStreamer::from_config(&config);
    let mut paired_streamer_raw = PairedRecordStreamer::with_pass_only(opt_c, config.vcf_pass_only);

    let mut table_buf = String::with_capacity(512);
    let mut vcf_scratch = String::with_capacity(512);
    let mut unimodal_used = false;
    let mut running_pos = i32::try_from(reg_start).unwrap_or(100);

    for row_idx in 0..num_rows {
        let pos_override = if clustered_positions {
            running_pos = running_pos.saturating_add(i32::from(cur.u8() % 16));
            Some(running_pos)
        } else {
            None
        };
        let allow_unimodal_here = enable_unimodal_once && row_idx == unimodal_target_row;
        let v1 = synth_variant(
            &mut cur,
            pos_override,
            allow_unimodal_here,
            &mut unimodal_used,
        );
        let sv_str = v1.sv.clone();

        // 1. Simple .var formatters (exact 36 and 38 columns).
        simple_36_into(&mut table_buf, &v1, &table_region, &table_sample, &sv_str);
        assert_tsv_cols("simple_36", &table_buf, 36);

        simple_38_into(&mut table_buf, &v1, &table_region, &table_sample, &sv_str);
        assert_tsv_cols("simple_38", &table_buf, 38);

        // 2. Amplicon .var formatters (exact 38 and 40 columns when debug == false,
        //    plus debug == true expansion).
        let gdamp = i32::from(cur.u8() % 8);
        let tlamp = gdamp.saturating_add(i32::from(cur.u8() % 8));
        let ncamp = i32::from(cur.u8() % 4);
        let ampflag = cur.bool();
        let av = AmpliconVariant {
            variant: v1.clone(),
            gdamp,
            tlamp,
            ncamp,
            ampflag,
            region: cur.utf8_capped(32, true),
            good_variants: if cur.bool() {
                vec![(Some(v1.clone()), "chr1:1-100".to_string())]
            } else {
                Vec::new()
            },
            bad_variants: if cur.bool() {
                vec![(None, "chr1:50-150".to_string())]
            } else {
                Vec::new()
            },
        };

        amplicon_38_into(&mut table_buf, &av, &table_region, &table_sample, false);
        assert_tsv_cols("amplicon_38", &table_buf, 38);

        amplicon_40_into(&mut table_buf, &av, &table_region, &table_sample, false);
        assert_tsv_cols("amplicon_40", &table_buf, 40);

        if row_idx == 0 {
            let extra_debug_cols = 1 + av.good_variants.len() + av.bad_variants.len();
            amplicon_38_into(&mut table_buf, &av, &table_region, &table_sample, true);
            assert_tsv_cols("amplicon_38_debug", &table_buf, 38 + extra_debug_cols);
            amplicon_40_into(&mut table_buf, &av, &table_region, &table_sample, true);
            assert_tsv_cols("amplicon_40_debug", &table_buf, 40 + extra_debug_cols);
        }

        // 3. Single / amplicon VCF event builder, direct formatter, and delayed streamer.
        let amp_opt = if cur.bool() {
            Some((gdamp, tlamp, ncamp, ampflag))
        } else {
            None
        };
        let custom_filters: Vec<String> = match cur.u8() % 3 {
            0 => Vec::new(),
            1 => vec!["q22.5".to_string(), "Bias".to_string()],
            _ => vec![{
                let f = cur.utf8_capped(16, false);
                if f.is_empty() {
                    "PASS".to_string()
                } else {
                    f
                }
            }],
        };
        let oddratio_tok = match cur.u8() % 4 {
            0 => "0",
            1 => "Inf",
            2 => "0.25",
            _ => "3.14159",
        };
        format_vcf_record_into(
            &mut vcf_scratch,
            &v1,
            &vcf_region,
            &config,
            &custom_filters,
            cur.metric_f64(1.0),
            oddratio_tok,
            amp_opt,
        );
        assert_vcf_data_record("single_vcf_direct", &vcf_scratch, 10, &config);

        // Fuzz SingleVcfRecord::from_parts_checked with arbitrary byte indices.
        let f_start = usize::from(cur.u16()) % (vcf_scratch.len().saturating_add(4));
        let f_end = usize::from(cur.u16()) % (vcf_scratch.len().saturating_add(4));
        if let Some(checked) = SingleVcfRecord::from_parts_checked(
            vcf_scratch.clone(),
            f_start,
            f_end,
            v1.start_pos,
            cur.bool(),
            cur.bool(),
            cur.bool(),
        ) {
            assert_vcf_data_record("single_vcf_checked", checked.as_str(), 10, &config);
            let _ = single_streamer_raw.push_event(SingleVcfEvent::Record(checked), |line| {
                assert_vcf_data_record("single_vcf_stream_checked", line, 10, &config);
                Ok(())
            });
        }

        let ev = build_single_vcf_event(&mut vcf_scratch, &v1, &vcf_region, &config, amp_opt);
        let _ = single_streamer_cfg.push_event(ev.clone(), |line| {
            assert_vcf_data_record("single_vcf_stream_cfg", line, 10, &config);
            Ok(())
        });
        let _ = single_streamer_raw.push_event(ev, |line| {
            assert_vcf_data_record("single_vcf_stream_raw", line, 10, &config);
            Ok(())
        });

        // 4. Somatic .var formatters (exact 55 and 61 columns) & Paired VCF streamer.
        //    Normal variant `v2` never enables unimodal (`false`) so `fisher2` and `fisher_total`
        //    always have support length <= 4097.
        let v2_opt = if cur.bool() {
            let mut v2 = synth_variant(&mut cur, Some(v1.start_pos), false, &mut unimodal_used);
            v2.ref_allele = v1.ref_allele.clone();
            v2.var_allele = v1.var_allele.clone();
            Some(v2)
        } else {
            None
        };
        let status = match cur.u8() % 9 {
            0 => "StrongSomatic".to_string(),
            1 => "LikelySomatic".to_string(),
            2 => "Germline".to_string(),
            3 => "LikelyLOH".to_string(),
            4 => "StrongLOH".to_string(),
            5 => "AFDiff".to_string(),
            6 => "SampleSpecific".to_string(),
            7 => "Deletion".to_string(),
            _ => cur.utf8_no_semi(16),
        };
        let sv2_str = match cur.u8() % 3 {
            0 => String::new(),
            1 => "0".to_string(),
            _ => cur.utf8_capped(16, false),
        };
        let end_is_variant2 = cur.bool();

        somatic_55_into(
            &mut table_buf,
            &v1,
            v2_opt.as_ref(),
            &table_region,
            &table_sample,
            &status,
            &sv_str,
            &sv2_str,
            end_is_variant2,
        );
        assert_tsv_cols("somatic_55", &table_buf, 55);

        somatic_61_into(
            &mut table_buf,
            &v1,
            v2_opt.as_ref(),
            &table_region,
            &table_sample,
            &status,
            &sv_str,
            &sv2_str,
            end_is_variant2,
        );
        assert_tsv_cols("somatic_61", &table_buf, 61);

        let som_var = SomaticVariant {
            variant1: v1,
            variant2: v2_opt,
            status,
            sv1: sv_str,
            sv2: sv2_str,
            end_is_variant2,
        };
        let paired_rec = format_paired_vcf_record(&som_var, &vcf_chr, &config);
        assert_vcf_data_record("paired_vcf_direct", paired_rec.as_str(), 11, &config);

        let _ = paired_streamer_cfg.push(paired_rec.clone(), |line| -> Result<(), ()> {
            assert_vcf_data_record("paired_vcf_stream_cfg", line, 11, &config);
            Ok(())
        });
        let _ = paired_streamer_raw.push(paired_rec, |line| -> Result<(), ()> {
            assert_vcf_data_record("paired_vcf_stream_raw", line, 11, &config);
            Ok(())
        });
    }

    let _ = single_streamer_cfg.finish(|line| {
        assert_vcf_data_record("single_vcf_finish_cfg", line, 10, &config);
        Ok(())
    });
    let _ = single_streamer_raw.finish(|line| {
        assert_vcf_data_record("single_vcf_finish_raw", line, 10, &config);
        Ok(())
    });
    let _ = paired_streamer_cfg.finish(|line| -> Result<(), ()> {
        assert_vcf_data_record("paired_vcf_finish_cfg", line, 11, &config);
        Ok(())
    });
    let _ = paired_streamer_raw.finish(|line| -> Result<(), ()> {
        assert_vcf_data_record("paired_vcf_finish_raw", line, 11, &config);
        Ok(())
    });

    // 5. Single & Paired VCF header writers (zero disk I/O via cleared reference path).
    let mut fuzz_chrom_lengths: HashMap<String, usize> = HashMap::with_capacity(4);
    let n_contigs = usize::from(cur.u8() % 5);
    for c_idx in 0..n_contigs {
        let cname = match cur.u8() % 6 {
            0 => "chr1".to_string(),
            1 => "chrX".to_string(),
            2 => "MT".to_string(),
            3 => format!("chr{}", (c_idx % 22) + 1),
            4 => vcf_chr.clone(),
            _ => cur.utf8_capped(24, true),
        };
        let clen = usize::from(cur.u16()).saturating_mul(64).saturating_add(1);
        fuzz_chrom_lengths.insert(cname, clen);
    }

    let mut hdr_cfg = config.clone();
    hdr_cfg.reference.clear();
    if (cur.u8() & 3) == 0 {
        hdr_cfg.sample = match cur.u8() % 4 {
            0 => "TUMOR A\u{0085}B|NORMAL C\u{2028}D\u{2029}E".to_string(),
            1 => "SAMPLE WITH SPACES\u{0085}\u{2028}\u{2029}".to_string(),
            2 => format!(
                "{}|{}",
                cur.utf8_capped(24, false),
                cur.utf8_capped(24, false)
            ),
            _ => cur.utf8_capped(48, true),
        };
    }

    let compat_off = hdr_cfg.compat.target == vardictrust_core::config::CompatTarget::Off;
    let mut hdr_buf: Vec<u8> = Vec::with_capacity(4096);
    if vardictrust_core::output::vcf::print_vcf_header(&mut hdr_buf, &hdr_cfg, &fuzz_chrom_lengths)
        .is_ok()
    {
        assert_vcf_header("single_vcf_header", &hdr_buf, 10, compat_off);
    }
    hdr_buf.clear();
    if vardictrust_core::output::vcf::print_paired_vcf_header(
        &mut hdr_buf,
        &hdr_cfg,
        &fuzz_chrom_lengths,
    )
    .is_ok()
    {
        assert_vcf_header("paired_vcf_header", &hdr_buf, 11, compat_off);
    }
});
