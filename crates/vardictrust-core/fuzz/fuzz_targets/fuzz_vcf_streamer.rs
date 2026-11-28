#![no_main]

use std::collections::HashMap;

use libfuzzer_sys::fuzz_target;
use vardictrust_core::config::{
    AmpliconConcordance, CompatProfile, Config, VcfThreshold, COMPAT_OFF, COMPAT_SWITCH_COUNT,
    COMPAT_VARDICTJAVA_1_8_3,
};
use vardictrust_core::output::table::{
    amplicon_row_into_compat, simple_row_into_compat, somatic_row_into_compat,
};
use vardictrust_core::output::vcf::{
    build_single_vcf_event, effective_vcf_cluster_dist, format_paired_vcf_record,
    format_vcf_record_into, PairedRecordStreamer, SingleRecordStreamer, SingleVcfEvent,
    SingleVcfRecord,
};
use vardictrust_core::postprocess::amplicon::AmpliconVariant;
use vardictrust_core::postprocess::simple::simple_post_process;
use vardictrust_core::postprocess::somatic::SomaticVariant;
use vardictrust_core::regions::Region;
use vardictrust_core::types::Pos1;
use vardictrust_core::variations::{Variant, Vars};

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
            9 => -((f64::from(self.u16()) / 65535.0) * max_norm.max(1.0)),
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
    /// IUPAC ambiguity codes, symbolic SV alleles (`<DUP>`, `<INV>`, `<DEL>`, `<INS>`, `<dup-N>`),
    /// and UTF-8 tokens.
    fn allele_256(&mut self) -> String {
        const BASES: &[u8; 5] = b"ACGTN";
        const IUPAC: &[u8; 10] = b"RYMKSWBDHV";
        match self.u8() % 13 {
            0 => String::new(),
            1 => "0".to_string(),
            2 => String::from(BASES[usize::from(self.u8()) % BASES.len()] as char),
            3 => "<DUP>".to_string(),
            4 => "<INV>".to_string(),
            5 => "<DEL>".to_string(),
            6 => "<INS>".to_string(),
            7 => format!("AC<dup-{}>GT", (self.u8() % 20) + 1),
            8 => String::from(IUPAC[usize::from(self.u8()) % IUPAC.len()] as char),
            9..=11 => {
                let len = ((usize::from(self.u8()) % MAX_ALLELE_BYTES) + 1).min(MAX_ALLELE_BYTES);
                let mut s = String::with_capacity(len);
                for _ in 0..len {
                    s.push(BASES[usize::from(self.u8()) % BASES.len()] as char);
                }
                s
            }
            _ => {
                let s = self.utf8_no_semi(MAX_ALLELE_BYTES).replace(',', "_");
                if s == "." {
                    "N".to_string()
                } else {
                    s
                }
            }
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

fn sanitize_allele_sw26(raw: &str) -> String {
    if raw.is_empty() || raw == "0" || raw == "." || is_symbolic_sv_alt(raw) {
        return raw.to_string();
    }
    let bytes = raw.as_bytes();
    let mut out = String::with_capacity(bytes.len());
    let mut i = 0usize;
    while i < bytes.len() {
        let b = bytes[i];
        if b == b'<' {
            if let Some(rel_close) = bytes[i + 1..].iter().position(|&c| c == b'>') {
                let inner = &raw[i + 1..i + 1 + rel_close];
                if !inner.is_empty()
                    && inner != "R,R"
                    && inner
                        .bytes()
                        .all(|c| c.is_ascii_alphanumeric() || c == b'-' || c == b'_')
                {
                    out.push_str(&raw[i..=i + 1 + rel_close]);
                    i += rel_close + 2;
                    continue;
                }
            }
        }
        let mapped = match b {
            b'A' | b'a' | b'W' | b'w' | b'M' | b'm' | b'R' | b'r' | b'D' | b'd' | b'H' | b'h'
            | b'V' | b'v' => 'A',
            b'C' | b'c' | b'S' | b's' | b'Y' | b'y' | b'B' | b'b' => 'C',
            b'G' | b'g' | b'K' | b'k' => 'G',
            b'T' | b't' | b'U' | b'u' => 'T',
            _ => 'N',
        };
        out.push(mapped);
        i += 1;
    }
    out
}

fn is_valid_sanitized_alt(alt: &str) -> bool {
    if alt == "." || is_symbolic_sv_alt(alt) {
        return true;
    }
    let bytes = alt.as_bytes();
    if bytes.is_empty() {
        return false;
    }
    let mut i = 0usize;
    while i < bytes.len() {
        let b = bytes[i];
        if b == b'<' {
            if let Some(rel_close) = bytes[i + 1..].iter().position(|&c| c == b'>') {
                let inner = &alt[i + 1..i + 1 + rel_close];
                if !inner.is_empty()
                    && inner != "R,R"
                    && inner
                        .bytes()
                        .all(|c| c.is_ascii_alphanumeric() || c == b'-' || c == b'_')
                {
                    i += rel_close + 2;
                    continue;
                }
            }
            return false;
        }
        if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N') {
            return false;
        }
        i += 1;
    }
    true
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

fn decode_sanitized_sample(sample_val: &str) -> String {
    if sample_val == "-" {
        return String::new();
    }
    let bytes = sample_val.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes.get(i) == Some(&b'%') {
            let Some(hex_slice) = bytes.get(i + 1..i + 3) else {
                panic!("truncated percent escape in INFO/SAMPLE: {sample_val:?}");
            };
            let hex_str = std::str::from_utf8(hex_slice)
                .unwrap_or_else(|e| panic!("non-ASCII percent escape in {sample_val:?}: {e}"));
            let byte = u8::from_str_radix(hex_str, 16)
                .unwrap_or_else(|e| panic!("invalid hex %{hex_str} in {sample_val:?}: {e}"));
            out.push(byte);
            i += 3;
        } else if let Some(&b) = bytes.get(i) {
            out.push(b);
            i += 1;
        } else {
            break;
        }
    }
    String::from_utf8(out)
        .unwrap_or_else(|e| panic!("decoded INFO/SAMPLE is not valid UTF-8 ({sample_val:?}): {e}"))
}

fn is_symbolic_sv_alt(alt: &str) -> bool {
    let Some(inner) = alt.strip_prefix('<').and_then(|s| s.strip_suffix('>')) else {
        return false;
    };
    !inner.is_empty()
        && inner != "R,R"
        && inner
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
}

fn render_expected_info_sample(raw_sample: &str, sanitize: bool) -> String {
    if !sanitize {
        return raw_sample
            .chars()
            .map(|c| if c.is_whitespace() { '_' } else { c })
            .collect();
    }
    if raw_sample.is_empty() {
        return "-".to_string();
    }
    if raw_sample == "-" {
        return "%2D".to_string();
    }
    if raw_sample == "." {
        return "%2E".to_string();
    }
    let mut buf = String::with_capacity(raw_sample.len());
    for c in raw_sample.chars() {
        match c {
            '%' => buf.push_str("%25"),
            ';' => buf.push_str("%3B"),
            '=' => buf.push_str("%3D"),
            ',' => buf.push_str("%2C"),
            c if c.is_whitespace() || c.is_control() => {
                let mut u8_buf = [0u8; 4];
                for &b in c.encode_utf8(&mut u8_buf).as_bytes() {
                    use std::fmt::Write as _;
                    let _ = write!(buf, "%{b:02X}");
                }
            }
            other => buf.push(other),
        }
    }
    buf
}

fn parse_nonneg_metric(kind: &str, label: &str, raw: &str, line: &str) {
    let val: f64 = raw
        .parse()
        .unwrap_or_else(|e| panic!("{kind} invalid float {label}={raw:?}: {e} in {line:?}"));
    assert!(
        val.is_finite() && val >= 0.0,
        "{kind} negative or non-finite {label}={raw:?} under nonnegative_subtracted_means: {line:?}"
    );
}

fn assert_vcf_data_record(kind: &str, line: &str, expected_cols: usize, config: &Config) {
    assert_tsv_cols(kind, line, expected_cols);
    if line.starts_with('#') {
        return;
    }
    let cols: Vec<&str> = line.split('\t').collect();
    let ref_col = cols.get(3).copied().unwrap_or("");
    let alt = cols.get(4).copied().unwrap_or("");
    let qual_str = cols.get(5).copied().unwrap_or("");
    let info = cols.get(7).copied().unwrap_or("");

    let expected_raw_sample = if expected_cols == 10 {
        config.sample.as_str()
    } else {
        config.sample.split('|').next().unwrap_or(&config.sample)
    };
    let rendered_sample =
        render_expected_info_sample(expected_raw_sample, config.compat.sanitize_sample_info);
    let sample_prefix = format!("SAMPLE={rendered_sample};");
    let info_after_sample = if expected_cols == 10 {
        info.strip_prefix(&sample_prefix)
            .unwrap_or_else(|| panic!("{kind} missing {sample_prefix:?} prefix in {info:?}"))
    } else {
        let (_, after_status) = info
            .split_once(';')
            .unwrap_or_else(|| panic!("{kind} missing STATUS; prefix in {info:?}"));
        after_status
            .strip_prefix(&sample_prefix)
            .unwrap_or_else(|| panic!("{kind} missing {sample_prefix:?} after STATUS in {info:?}"))
    };

    // Switch 26 (`sanitize_iupac_alt_alleles`): no IUPAC ambiguity codes outside `ACGTN`/`.`/`<...>`.
    if config.compat.sanitize_iupac_alt_alleles {
        assert!(
            is_valid_sanitized_alt(alt),
            "{kind} leaked unsanitized ALT={alt:?} under sanitize_iupac_alt_alleles: {line:?}"
        );
    }

    // Switch 23 (`nonnegative_subtracted_means`): non-negative VCF QUAL, MQ, NM.
    if config.compat.nonnegative_subtracted_means {
        let qual: i32 = qual_str
            .parse()
            .unwrap_or_else(|e| panic!("{kind} invalid integer QUAL {qual_str:?}: {e}"));
        assert!(
            qual >= 0,
            "{kind} negative QUAL={qual} under nonnegative_subtracted_means: {line:?}"
        );
        if expected_cols == 10 {
            for kv in info_after_sample.split(';') {
                if let Some(v) = kv.strip_prefix("QUAL=") {
                    parse_nonneg_metric(kind, "INFO/QUAL", v, line);
                } else if let Some(v) = kv.strip_prefix("MQ=") {
                    parse_nonneg_metric(kind, "INFO/MQ", v, line);
                } else if let Some(v) = kv.strip_prefix("NM=") {
                    parse_nonneg_metric(kind, "INFO/NM", v, line);
                }
            }
        } else {
            let fmt_fields: Vec<&str> = cols.get(8).copied().unwrap_or("").split(':').collect();
            let qual_idx = fmt_fields.iter().position(|&k| k == "QUAL");
            let mq_idx = fmt_fields.iter().position(|&k| k == "MQ");
            let nm_idx = fmt_fields.iter().position(|&k| k == "NM");
            for &s_col in &[9usize, 10usize] {
                let sample_str = cols.get(s_col).copied().unwrap_or("");
                let s_fields: Vec<&str> = sample_str.split(':').collect();
                assert!(
                    s_fields.len() >= fmt_fields.len(),
                    "{kind} sample col {s_col} has fewer fields ({}) than FORMAT ({}): {line:?}",
                    s_fields.len(),
                    fmt_fields.len()
                );
                let shift = s_fields.len() - fmt_fields.len();
                if let Some(q_v) = qual_idx.and_then(|i| s_fields.get(i + shift).copied()) {
                    parse_nonneg_metric(kind, "FORMAT/QUAL", q_v, line);
                }
                if let Some(mq_v) = mq_idx.and_then(|i| s_fields.get(i + shift).copied()) {
                    parse_nonneg_metric(kind, "FORMAT/MQ", mq_v, line);
                }
                if let Some(nm_v) = nm_idx.and_then(|i| s_fields.get(i + shift).copied()) {
                    parse_nonneg_metric(kind, "FORMAT/NM", nm_v, line);
                }
            }
        }
    }

    // Switch 30 (`strict_vcf_ad_number_r`): `len(AD) == 1 + n_alt` (`1` on `.`, `2` on single concrete/SV ALT).
    if config.compat.strict_vcf_ad_number_r {
        let n_alt = if alt == "." {
            0usize
        } else {
            alt.split(',').count()
        };
        let expected_ad_len = 1 + n_alt;
        let ad_idx = if expected_cols == 10 { 3usize } else { 5usize };
        let sample_indices: &[usize] = if expected_cols == 10 { &[9] } else { &[9, 10] };
        for &s_col in sample_indices {
            let sample_str = cols.get(s_col).copied().unwrap_or("");
            let ad_str = sample_str.split(':').nth(ad_idx).unwrap_or("");
            let ad_tokens: Vec<&str> = ad_str.split(',').collect();
            assert_eq!(
                ad_tokens.len(),
                expected_ad_len,
                "{kind} FORMAT/AD token count mismatch at col {s_col} (alt={alt:?}): got {ad_str:?}, expected {expected_ad_len}"
            );
            for tok in ad_tokens {
                assert!(
                    tok.parse::<i64>().is_ok(),
                    "{kind} non-integer AD token {tok:?} in {ad_str:?}"
                );
            }
        }
    }

    // Switch 31 (`omit_nonsv_vcf_end`): `END` iff symbolic SV unless `--no-end` (`vcf_no_end_explicit`).
    if config.compat.omit_nonsv_vcf_end {
        let has_end = info_after_sample
            .split(';')
            .any(|kv| kv.starts_with("END="));
        if expected_cols == 10 || config.compat.canonical_sv_vcf_anchoring {
            let want_end = is_symbolic_sv_alt(alt) && !config.vcf_no_end_explicit;
            assert_eq!(
                has_end, want_end,
                "{kind} INFO/END mismatch under omit_nonsv_vcf_end (alt={alt:?}, no_end={}): {info:?}",
                config.vcf_no_end_explicit
            );
        } else if !is_symbolic_sv_alt(alt) || config.vcf_no_end_explicit {
            assert!(
                !has_end,
                "{kind} emitted non-SV or --no-end INFO/END under omit_nonsv_vcf_end (alt={alt:?}): {info:?}"
            );
        }
    }

    // Switch 38 (`indel_only_shift3`): `SHIFT3 == 0` whenever `|REF| == |ALT|` (non-dot, non-symbolic).
    if config.compat.indel_only_shift3
        && !ref_col.is_empty()
        && alt != "."
        && !alt.contains('<')
        && ref_col.len() == alt.len()
    {
        let shift3_val: i32 = info_after_sample
            .split(';')
            .find_map(|kv| kv.strip_prefix("SHIFT3="))
            .and_then(|v| v.parse().ok())
            .unwrap_or_else(|| panic!("{kind} missing integer SHIFT3= in {info:?}"));
        assert_eq!(
            shift3_val, 0,
            "{kind} non-zero SHIFT3={shift3_val} on equal-length allele ref={ref_col:?} alt={alt:?}: {line:?}"
        );
    }

    // Switch 32 (`canonical_vcf_gt`): sorted unphased `GT` (`0/0`, `0/1`, `1/1`, never `1/0`).
    if config.compat.canonical_vcf_gt {
        let sample_indices: &[usize] = if expected_cols == 10 { &[9] } else { &[9, 10] };
        let mut any_vd_pos = false;
        let mut any_gt_one = false;
        for &s_col in sample_indices {
            let sample_str = cols.get(s_col).copied().unwrap_or("");
            let mut s_parts = sample_str.split(':');
            let gt = s_parts.next().unwrap_or("");
            let _dp = s_parts.next();
            let vd: i32 = s_parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            assert_ne!(
                gt, "1/0",
                "{kind} emitted unphased 1/0 under canonical_vcf_gt: {sample_str:?}"
            );
            assert!(
                matches!(gt, "0/0" | "0/1" | "1/1"),
                "{kind} unexpected GT {gt:?}: {sample_str:?}"
            );
            if alt == "." {
                assert_eq!(
                    gt, "0/0",
                    "{kind} dot-ALT must have GT=0/0 under canonical_vcf_gt: {sample_str:?}"
                );
            }
            if vd > 0 {
                any_vd_pos = true;
            }
            if gt.contains('1') {
                any_gt_one = true;
            }
        }
        if alt != "." && any_vd_pos {
            assert!(
                any_gt_one,
                "{kind} concrete ALT={alt:?} with VD>0 has all 0/0 GTs under canonical_vcf_gt: {line:?}"
            );
        }
    }

    if !config.compat.sanitize_sample_info {
        return;
    }
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
            && !sample_val
                .chars()
                .any(|c| c.is_whitespace() || c.is_control()),
        "{kind} has unescaped VCF delimiter, whitespace, or control char in INFO/SAMPLE: {sample_val:?}"
    );
    assert_eq!(
        decode_sanitized_sample(sample_val),
        expected_raw_sample,
        "{kind} INFO/SAMPLE round-trip mismatch: encoded={sample_val:?}, expected_raw={expected_raw_sample:?}"
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
        let emit_end = if config.compat.omit_nonsv_vcf_end {
            is_symbolic_sv_alt(alt) && !config.vcf_no_end_explicit
        } else if config.compat.canonical_sv_vcf_anchoring && is_symbolic_sv_alt(alt) {
            !config.vcf_no_end_explicit
        } else {
            !config.omit_end_tag
        };
        let base_end = 26 + usize::from(emit_end);
        if entries.iter().any(|kv| kv.starts_with("GDAMP=")) {
            let amp_sv_extra =
                usize::from(config.compat.canonical_sv_vcf_anchoring && is_symbolic_sv_alt(alt))
                    * 2;
            base_end + 4 + amp_sv_extra
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
        let amp_extra = if entries.iter().any(|kv| kv.starts_with("GDAMP=")) {
            4
        } else {
            0
        };
        let sv_extra = if entries.iter().any(|kv| kv.starts_with("SVTYPE=")) {
            2 + usize::from(!config.vcf_no_end_explicit)
        } else {
            0
        };
        let flag_extra = usize::from(
            config.vcf_paired_strict
                && matches!(entries.last().copied(), Some("P0.01Likely" | "InDelLikely")),
        );
        13 + amp_extra + sv_extra + flag_extra
    };
    assert_eq!(
        entries.len(),
        expected_keys,
        "{kind} INFO key count mismatch: got {}, expected {expected_keys} in {info:?}",
        entries.len()
    );
}

fn assert_vcf_header(
    kind: &str,
    hdr_bytes: &[u8],
    expected_chrom_cols: usize,
    compat_non_default: bool,
    preset: Option<vardictrust_core::config::WorkflowPreset>,
) {
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
    // S-08 (PR #168): `##VarDictRustCompat=off` is emitted iff `compat != COMPAT_VARDICTJAVA_1_8_3`.
    if compat_non_default {
        assert_eq!(
            lines
                .iter()
                .filter(|&&l| l == "##VarDictRustCompat=off")
                .count(),
            1,
            "{kind} must emit ##VarDictRustCompat=off once when compat != COMPAT_VARDICTJAVA_1_8_3"
        );
        assert_eq!(lines.get(2), Some(&"##VarDictRustCompat=off"));
        if let Some(p) = preset {
            let expected_preset_line = format!("##VarDictRustPreset={p}");
            assert_eq!(
                lines
                    .iter()
                    .filter(|&&l| l == expected_preset_line.as_str())
                    .count(),
                1,
                "{kind} must emit {expected_preset_line} once when preset={p}"
            );
            assert_eq!(lines.get(3), Some(&expected_preset_line.as_str()));
        } else {
            assert!(
                !lines.iter().any(|l| l.starts_with("##VarDictRustPreset=")),
                "{kind} must not emit ##VarDictRustPreset= when preset is None"
            );
        }
    } else {
        assert!(
            !lines.iter().any(|l| l.starts_with("##VarDictRustCompat")),
            "{kind} must not emit ##VarDictRustCompat under COMPAT_VARDICTJAVA_1_8_3"
        );
        assert!(
            !lines.iter().any(|l| l.starts_with("##VarDictRustPreset=")),
            "{kind} must not emit ##VarDictRustPreset= under COMPAT_VARDICTJAVA_1_8_3"
        );
    }
}

fn synth_variant(
    cur: &mut ByteCursor<'_>,
    pos_override: Option<i32>,
    allow_unimodal: bool,
    unimodal_used: &mut bool,
    compat: &CompatProfile,
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
    if compat.sanitize_iupac_alt_alleles {
        v.ref_allele = sanitize_allele_sw26(&v.ref_allele);
        v.var_allele = sanitize_allele_sw26(&v.var_allele);
        if v.var_allele.is_empty() || v.var_allele == "0" {
            v.var_allele = "N".to_string();
        }
    }
    debug_assert!(v.ref_allele.len() <= 256 && v.var_allele.len() <= 256);

    v.bias = match cur.u8() % 7 {
        0 => "0".to_string(),
        1 => "0;0".to_string(),
        2 => "2;1".to_string(),
        3 => "2;2".to_string(),
        4 => "1;2".to_string(),
        5 => "0;2".to_string(),
        _ => cur.utf8_capped(12, false).replace(':', "_"),
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

    // Mirror upstream `ToVarsBuilder` / `Variant::normalize_postprocess` invariants on direct
    // VCF formatter inputs when the corresponding switches are active (while Part 6 below
    // exercises `ToVarsBuilder` -> `simple_post_process` -> `SingleRecordStreamer` on raw inputs).
    if compat.nonnegative_subtracted_means {
        let clamp_nonneg = |x: f64| if x.is_finite() && x > 0.0 { x } else { 0.0 };
        v.mean_position = clamp_nonneg(v.mean_position);
        v.mean_quality = clamp_nonneg(v.mean_quality);
        v.mean_mapping_quality = clamp_nonneg(v.mean_mapping_quality);
        v.nm = clamp_nonneg(v.nm);
    }
    if compat.indel_only_shift3
        && !v.ref_allele.is_empty()
        && !v.var_allele.contains('<')
        && v.ref_allele.len() == v.var_allele.len()
    {
        v.shift3 = 0;
    }
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

    // Exercise pure `COMPAT_VARDICTJAVA_1_8_3`, pure `COMPAT_OFF`, single-switch, leave-one-out,
    // and full 38-bit `CompatProfile::from_mask(mask64)` across all 38 switches (`0..38`).
    let fuzz_compat = {
        let sel = cur.u8();
        let mask64 = cur.u64();
        match sel % 6 {
            0 => COMPAT_VARDICTJAVA_1_8_3,
            1 => COMPAT_OFF,
            2 => CompatProfile::only_switch((mask64 as usize) % COMPAT_SWITCH_COUNT)
                .unwrap_or(COMPAT_OFF),
            3 => CompatProfile::without_switch((mask64 as usize) % COMPAT_SWITCH_COUNT)
                .unwrap_or(COMPAT_OFF),
            _ => CompatProfile::from_mask(mask64),
        }
    };
    let fuzz_preset = if fuzz_compat != COMPAT_VARDICTJAVA_1_8_3 {
        use vardictrust_core::config::WorkflowPreset;
        match cur.u8() % 5 {
            0 => None,
            1 => Some(WorkflowPreset::Wgs),
            2 => Some(WorkflowPreset::Wes),
            3 => Some(WorkflowPreset::Amplicon),
            _ => Some(WorkflowPreset::Ctdna),
        }
    } else {
        None
    };

    let config = Config {
        fisher: (mode_bits & 2) != 0,
        vcf_pass_only: (mode_bits & 4) != 0,
        vcf_paired_strict: (mode_bits & 8) != 0,
        omit_end_tag: (mode_bits & 16) != 0,
        vcf_no_end_explicit: cur.bool(),
        vcf_pstd_filter: (mode_bits & 32) != 0,
        pileup: (mode_bits & 64) != 0,
        af_threshold_truthy_zero: cur.bool(),
        af_threshold: cur.metric_f64(0.5),
        af_threshold_raw: match cur.u8() % 6 {
            0 => None,
            1 => Some("0".to_string()),
            2 => Some("0.0".to_string()),
            3 => Some("0.050".to_string()),
            4 => Some("1e-5".to_string()),
            _ => Some(cur.utf8_capped(48, true)),
        },
        vcf_min_splitreads: i32::from(cur.u8() % 8) - 1,
        vcf_min_splitreads_explicit: cur.bool(),
        crispr_cutting_site: match cur.u8() % 3 {
            0 => None,
            1 => Some(0),
            _ => Some(cur.i32()),
        },
        sample: match cur.u8() % 12 {
            0 => "SAMPLE".to_string(),
            1 => "TUMOR|NORMAL".to_string(),
            2 => "TUMOR|".to_string(),
            3 => String::new(),
            4 => ".".to_string(),
            5 => "..|.".to_string(),
            6 => "a;FOO=1,bar%20baz\tq|NORM;X=2".to_string(),
            7 => "-".to_string(),
            8 => "%2D|_".to_string(),
            9 => "a b;T=1|a_b;T=1".to_string(),
            10 => "\0".to_string(),
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
        compat: fuzz_compat,
        preset: fuzz_preset,
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
            &config.compat,
        );
        let sv_str = v1.sv.clone();

        let exact_tie = config.compat.exact_binary_half_even;
        // 1. Simple .var formatters (exact 36 and 38 columns).
        simple_row_into_compat(
            &mut table_buf,
            &v1,
            &table_region,
            &table_sample,
            &sv_str,
            false,
            exact_tie,
        );
        assert_tsv_cols("simple_36", &table_buf, 36);

        simple_row_into_compat(
            &mut table_buf,
            &v1,
            &table_region,
            &table_sample,
            &sv_str,
            true,
            exact_tie,
        );
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

        amplicon_row_into_compat(
            &mut table_buf,
            &av,
            &table_region,
            &table_sample,
            false,
            false,
            exact_tie,
        );
        assert_tsv_cols("amplicon_38", &table_buf, 38);

        amplicon_row_into_compat(
            &mut table_buf,
            &av,
            &table_region,
            &table_sample,
            false,
            true,
            exact_tie,
        );
        assert_tsv_cols("amplicon_40", &table_buf, 40);

        if row_idx == 0 {
            let extra_debug_cols = 1 + av.good_variants.len() + av.bad_variants.len();
            amplicon_row_into_compat(
                &mut table_buf,
                &av,
                &table_region,
                &table_sample,
                true,
                false,
                exact_tie,
            );
            assert_tsv_cols("amplicon_38_debug", &table_buf, 38 + extra_debug_cols);
            amplicon_row_into_compat(
                &mut table_buf,
                &av,
                &table_region,
                &table_sample,
                true,
                true,
                exact_tie,
            );
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
        let v2_opt = if cur.bool() {
            let mut v2 = synth_variant(
                &mut cur,
                Some(v1.start_pos),
                false,
                &mut unimodal_used,
                &config.compat,
            );
            v2.ref_allele = v1.ref_allele.clone();
            v2.var_allele = v1.var_allele.clone();
            if config.compat.indel_only_shift3
                && !v2.ref_allele.is_empty()
                && !v2.var_allele.contains('<')
                && v2.ref_allele.len() == v2.var_allele.len()
            {
                v2.shift3 = 0;
            }
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

        somatic_row_into_compat(
            &mut table_buf,
            &v1,
            v2_opt.as_ref(),
            &table_region,
            &table_sample,
            &status,
            &sv_str,
            &sv2_str,
            end_is_variant2,
            false,
            exact_tie,
        );
        assert_tsv_cols("somatic_55", &table_buf, 55);

        somatic_row_into_compat(
            &mut table_buf,
            &v1,
            v2_opt.as_ref(),
            &table_region,
            &table_sample,
            &status,
            &sv_str,
            &sv2_str,
            end_is_variant2,
            true,
            exact_tie,
        );
        assert_tsv_cols("somatic_61", &table_buf, 61);

        let synth_amp_conc = |c: &mut ByteCursor<'_>| -> Option<AmpliconConcordance> {
            if c.bool() {
                let gd = i32::from(c.u8() % 6);
                let tl = gd.saturating_add(i32::from(c.u8() % 6));
                let nc = i32::from(c.u8() % 4);
                Some(AmpliconConcordance {
                    gdamp: gd,
                    tlamp: tl,
                    ncamp: nc,
                    ampflag: c.bool(),
                })
            } else {
                None
            }
        };
        let tumor_amplicon = synth_amp_conc(&mut cur);
        let normal_amplicon = synth_amp_conc(&mut cur);
        let som_var = SomaticVariant {
            variant1: v1,
            variant2: v2_opt,
            status,
            sv1: sv_str,
            sv2: sv2_str,
            end_is_variant2,
            tumor_amplicon,
            normal_amplicon,
            v2_ref_variant: None,
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

    // S-08: `##VarDictRustCompat=off` is emitted iff `compat != COMPAT_VARDICTJAVA_1_8_3`.
    let compat_non_default = hdr_cfg.compat != COMPAT_VARDICTJAVA_1_8_3;
    let mut hdr_buf: Vec<u8> = Vec::with_capacity(4096);
    if vardictrust_core::output::vcf::print_vcf_header(&mut hdr_buf, &hdr_cfg, &fuzz_chrom_lengths)
        .is_ok()
    {
        assert_vcf_header(
            "single_vcf_header",
            &hdr_buf,
            10,
            compat_non_default,
            hdr_cfg.preset,
        );
    }
    hdr_buf.clear();
    if vardictrust_core::output::vcf::print_paired_vcf_header(
        &mut hdr_buf,
        &hdr_cfg,
        &fuzz_chrom_lengths,
    )
    .is_ok()
    {
        assert_vcf_header(
            "paired_vcf_header",
            &hdr_buf,
            11,
            compat_non_default,
            hdr_cfg.preset,
        );
    }

    // 6. End-to-end `ToVarsBuilder` -> `simple_post_process` -> `SingleRecordStreamer` +
    //    `simple_row_into_compat` invariants across all 38 `CompatProfile` switches (injecting
    //    real IUPAC ambiguity alleles, homopolymer runs for `SHIFT3`, and negative means).
    {
        let span_u32 = reg_end.saturating_sub(reg_start).clamp(16, 512);
        let eff_end = reg_start.saturating_add(span_u32);
        let w_region = Region {
            chr: "chr1".to_string(),
            start: Pos1(reg_start),
            end: Pos1(eff_end),
            gene: "W_FUZZ".to_string(),
            amplicon: None,
        };
        let ref_len = (span_u32 as usize).saturating_add(1);
        let ref_seq: Vec<u8> = (0..ref_len)
            .map(|i| {
                if i % 16 < 6 {
                    b'A' // homopolymer run to trigger non-zero MSI / raw shift3
                } else {
                    b"ACGT"[(usize::from(cur.u8()) + i) % 4]
                }
            })
            .collect();
        let mut reference = vardictrust_core::reference::Reference::new();
        reference
            .loaded_regions
            .push(vardictrust_core::reference::LoadedRegion {
                chr: "chr1".to_string(),
                sequence_start: reg_start,
                sequence_end: eff_end,
            });
        reference.region_sequences.push(ref_seq.clone());
        let mut w_chr_lengths = HashMap::with_capacity(1);
        w_chr_lengths.insert("chr1".to_string(), eff_end as usize);

        let mut vd = vardictrust_core::variations::VariationData::new();
        let n_sites = ((usize::from(cur.u8()) % 8) + 2).min(8);
        const ALT_POOL: &[char; 8] = &['A', 'C', 'G', 'T', 'R', 'Y', 'W', 'S'];
        for _ in 0..n_sites {
            let off = u32::from(cur.u16()) % (span_u32.max(1));
            let pos_i32 = i32::try_from(reg_start.saturating_add(off)).unwrap_or(1);
            let raw_dp = i32::from(cur.u8() % 32);
            vd.ref_coverage.insert(pos_i32, raw_dp);

            let ref_ch = ref_seq.get(off as usize).copied().unwrap_or(b'A') as char;
            let alt_ch = {
                let cand = ALT_POOL[usize::from(cur.u8()) % ALT_POOL.len()];
                if cand == ref_ch {
                    if ref_ch == 'A' {
                        'G'
                    } else {
                        'A'
                    }
                } else {
                    cand
                }
            };
            let mut vmap = vardictrust_core::variations::VariationMap::new();
            let mut ref_var = vardictrust_core::variations::Variation::new();
            let rfwd = i32::from(cur.u8() % 16);
            let rrev = i32::from(cur.u8() % 16);
            ref_var.cnt = rfwd + rrev;
            ref_var.fwd = rfwd;
            ref_var.rev = rrev;
            ref_var.hicnt = i32::from(cur.u8() % 20) - 2;
            ref_var.locnt = i32::from(cur.u8() % 4);
            ref_var.extracnt = i32::from(cur.u8() % 20);
            ref_var.pmean = f64::from(ref_var.cnt.max(1) * 15);
            ref_var.qmean = f64::from(ref_var.cnt.max(1) * 35);
            ref_var.mqmean = f64::from(ref_var.cnt.max(1) * 60);
            vmap.map.insert(ref_ch.to_string(), ref_var);

            let mut alt_var = vardictrust_core::variations::Variation::new();
            let afwd = i32::from(cur.u8() % 20) + 2;
            let arev = i32::from(cur.u8() % 20) + 1;
            alt_var.cnt = afwd + arev;
            alt_var.fwd = afwd;
            alt_var.rev = arev;
            alt_var.hicnt = alt_var.cnt;
            alt_var.extracnt = i32::from(cur.u8() % (alt_var.cnt as u8).max(1));
            if (cur.u8() & 7) == 0 {
                alt_var.pmean = -5.0;
                alt_var.qmean = -10.0;
                alt_var.mqmean = -3.0;
                alt_var.nm = -2.0;
            } else {
                alt_var.pmean = f64::from(alt_var.cnt * 20);
                alt_var.qmean = f64::from(alt_var.cnt * 36);
                alt_var.mqmean = f64::from(alt_var.cnt * 60);
            }
            alt_var.pstd = true;
            alt_var.qstd = true;
            vmap.map.insert(alt_ch.to_string(), alt_var);
            vd.non_insertion_variants.insert(pos_i32, vmap);
        }

        let mut w_cfg = config.clone();
        w_cfg.emit_intermediate = true;
        w_cfg.pileup = cur.bool();
        if w_cfg.compat.coordinate_sorted_vars {
            w_cfg.compat.suppress_pileup_sentinels = true;
        }
        let builder = vardictrust_core::tovars::ToVarsBuilder::new(
            &reference,
            &w_region,
            &w_chr_lengths,
            &w_cfg,
        );
        let aligned: vardictrust_core::javamap::JavaHashMap<i32, Vars> =
            builder.build_aligned_variants(&vd);

        let check_reconciled = |v: &Variant, is_ref: bool| {
            assert!(
                0 <= v.pos_coverage && v.pos_coverage <= v.total_pos_cov,
                "VD out of [0, DP]: vd={}, dp={}, is_ref={is_ref}",
                v.pos_coverage,
                v.total_pos_cov
            );
            if is_ref {
                assert!(
                    0 <= v.ref_fwd_cov
                        && 0 <= v.ref_rev_cov
                        && v.ref_fwd_cov.saturating_add(v.ref_rev_cov) <= v.total_pos_cov,
                    "ref_variant RD out of [0, DP]: rf={}, rr={}, dp={}",
                    v.ref_fwd_cov,
                    v.ref_rev_cov,
                    v.total_pos_cov
                );
            }
            assert!(0.0 <= v.frequency && v.frequency <= 1.0);
            assert!(0.0 <= v.extra_freq && v.extra_freq <= 1.0);
            assert!(0 <= v.hicnt && v.hicnt <= v.hicov);
            assert!(0.0 <= v.hi_freq && v.hi_freq <= 1.0);
        };

        if w_cfg.compat.reconcile_vd_dp || w_cfg.compat.reconcile_rd_vd_dp {
            for (_, vars) in aligned.iter() {
                for v in &vars.variants {
                    check_reconciled(v, false);
                }
                if let Some(ref rv) = vars.ref_variant {
                    check_reconciled(rv, true);
                }
            }
        }

        let emitted = simple_post_process(&w_cfg, &w_region, &aligned, &vd.splice);
        let mut prev_pos = i32::MIN;
        let mut e2e_streamer = SingleRecordStreamer::from_config(&w_cfg);
        for v in &emitted {
            if w_cfg.compat.reconcile_vd_dp || w_cfg.compat.reconcile_rd_vd_dp {
                check_reconciled(v, v.var_type.is_empty());
            }
            simple_row_into_compat(
                &mut table_buf,
                v,
                &w_region,
                "W_SMPL",
                &v.sv,
                false,
                w_cfg.compat.exact_binary_half_even,
            );
            assert_tsv_cols("w_simple_36", &table_buf, 36);
            if w_cfg.compat.coordinate_sorted_vars {
                let col_pos: i32 = table_buf
                    .split('\t')
                    .nth(3)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(v.start_pos);
                assert!(
                    col_pos >= prev_pos,
                    "coordinate_sorted_vars .var Start regressed: {col_pos} < {prev_pos}"
                );
                prev_pos = col_pos;
            }
            let e2e_ev = build_single_vcf_event(&mut vcf_scratch, v, &w_region, &w_cfg, None);
            let _ = e2e_streamer.push_event(e2e_ev, |line| {
                assert_vcf_data_record("e2e_single_vcf", line, 10, &w_cfg);
                Ok(())
            });
        }
        let _ = e2e_streamer.finish(|line| {
            assert_vcf_data_record("e2e_single_vcf_finish", line, 10, &w_cfg);
            Ok(())
        });
    }
});
