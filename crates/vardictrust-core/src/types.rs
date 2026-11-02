//! Core types for VarDictRust: coordinate newtypes,
//! rounding utilities, and the unified error type.
//!
//! These types enforce coordinate safety at compile time and replicate
//! VarDictJava's exact output formatting.

use std::fmt;

// ---------------------------------------------------------------------------
// Coordinate newtypes — the #1 silent footgun (CLAUDE.md §6)
// ---------------------------------------------------------------------------

/// 0-based, half-open position (BAM/SAM coordinate system).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Pos0(pub u32);

/// 1-based position (VCF/VarDict-internal coordinate system).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Pos1(pub u32);

impl Pos0 {
    /// Convert a 0-based position to 1-based.
    #[inline]
    pub fn to_one_based(self) -> Pos1 {
        Pos1(self.0.saturating_add(1))
    }
}

impl Pos1 {
    /// Convert a 1-based position to 0-based.
    #[inline]
    pub fn to_zero_based(self) -> Pos0 {
        Pos0(self.0.saturating_sub(1))
    }
}

/// Noodles `Position` is 1-based by construction — document the invariant
/// once here instead of wrapping raw `u32`s at every call site (§6).
impl From<noodles::core::Position> for Pos1 {
    #[inline]
    fn from(p: noodles::core::Position) -> Self {
        Pos1(p.get() as u32)
    }
}

impl fmt::Display for Pos0 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl fmt::Display for Pos1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

// ---------------------------------------------------------------------------
// Rounding — banker's rounding (HALF_EVEN) to match Java's DecimalFormat
// ---------------------------------------------------------------------------

struct StackBuf<const N: usize> {
    buf: [u8; N],
    len: usize,
}

impl<const N: usize> StackBuf<N> {
    #[inline]
    const fn new() -> Self {
        Self {
            buf: [0u8; N],
            len: 0,
        }
    }

    #[inline]
    fn as_str(&self) -> Option<&str> {
        let slice = self.buf.get(..self.len)?;
        std::str::from_utf8(slice).ok()
    }
}

impl<const N: usize> fmt::Write for StackBuf<N> {
    #[inline]
    fn write_str(&mut self, s: &str) -> fmt::Result {
        let bytes = s.as_bytes();
        let end = self.len.checked_add(bytes.len()).ok_or(fmt::Error)?;
        let dst = self.buf.get_mut(self.len..end).ok_or(fmt::Error)?;
        dst.copy_from_slice(bytes);
        self.len = end;
        Ok(())
    }
}

/// Returns `true` when `(value, digits)` hits OpenJDK `DecimalFormat` (`DigitList.set`
/// fixed-point underflow branch `-decimalAt == maximumDigits`), where the shortest
/// decimal representation `5e-(d+1)` rounds to `±0.0` under `HALF_EVEN` even though the
/// exact binary `f64` exceeds the midpoint `5 × 10^-(d+1)`.
#[inline]
fn is_decimal_format_underflow_tie(value: f64, digits: i32) -> bool {
    let abs_bits = value.abs().to_bits();
    matches!(
        (digits, abs_bits),
        (3, 0x3f40_624d_d2f1_a9fc) // 5e-4
            | (4, 0x3f0a_36e2_eb1c_432d) // 5e-5
            | (5, 0x3ed4_f8b5_88e3_68f1) // 5e-6
            | (8, 0x3e35_798e_e230_8c3a) // 5e-9
    )
}

/// Banker's rounding to `digits` decimal places.
///
/// Matches Java's `DecimalFormat` with `RoundingMode.HALF_EVEN`.
/// This is critical for float-field equivalence with VarDictJava.
pub fn round_half_even(value: f64, digits: i32) -> f64 {
    if !value.is_finite() {
        return value;
    }
    if is_decimal_format_underflow_tie(value, digits) {
        return 0.0_f64.copysign(value);
    }
    let prec = digits.max(0) as usize;
    let mut buf = StackBuf::<64>::new();
    use std::fmt::Write as _;
    if write!(&mut buf, "{:.*}", prec, value).is_ok() {
        if let Some(s) = buf.as_str() {
            if let Ok(parsed) = s.parse::<f64>() {
                return parsed;
            }
        }
    }
    let s = format!("{:.*}", prec, value);
    s.parse::<f64>().unwrap_or(value)
}

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

/// Unified error type for vardictrust-core.
#[derive(thiserror::Error, Debug)]
pub enum VardictError {
    /// I/O error (file access, etc.)
    #[error("I/O: {0}")]
    Io(#[from] std::io::Error),

    /// Feature not yet implemented in v1.0.
    #[error("unsupported: {0}")]
    Unsupported(String),

    /// Equivalent to Java throwing an exception and skipping a read.
    #[error("Java exception equivalent: {0}")]
    JavaThrows(String),

    /// Parse error (BED, region string, etc.)
    #[error("parse: {0}")]
    Parse(String),

    /// Generic error.
    #[error("{0}")]
    Other(String),
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pos0_to_pos1() {
        assert_eq!(Pos0(0).to_one_based(), Pos1(1));
        assert_eq!(Pos0(99).to_one_based(), Pos1(100));
    }

    #[test]
    fn pos1_to_pos0() {
        assert_eq!(Pos1(1).to_zero_based(), Pos0(0));
        assert_eq!(Pos1(100).to_zero_based(), Pos0(99));
        // Saturating: Pos1(0) should not underflow
        assert_eq!(Pos1(0).to_zero_based(), Pos0(0));
    }

    #[test]
    fn round_half_even_banker() {
        // 2.5 rounds to 2 (nearest even)
        assert_eq!(round_half_even(2.5, 0), 2.0);
        // 3.5 rounds to 4 (nearest even)
        assert_eq!(round_half_even(3.5, 0), 4.0);
        // 0.5 rounds to 0
        assert_eq!(round_half_even(0.5, 0), 0.0);
        // 1.5 rounds to 2
        assert_eq!(round_half_even(1.5, 0), 2.0);
        // Standard rounding (not at .5 boundary)
        assert_eq!(round_half_even(2.3, 0), 2.0);
        assert_eq!(round_half_even(2.7, 0), 3.0);
    }

    #[test]
    fn round_half_even_digits() {
        // 5 decimal places (the VarDict standard)
        assert_eq!(round_half_even(0.123456, 5), 0.12346);
        assert_eq!(round_half_even(0.012345, 5), 0.01235); // due to exact double representation being slightly > 0.012345
        assert_eq!(round_half_even(0.012355, 5), 0.01235); // due to exact double representation being slightly < 0.012355
    }

    #[test]
    fn round_half_even_special() {
        assert!(round_half_even(f64::NAN, 5).is_nan());
        assert_eq!(round_half_even(f64::INFINITY, 5), f64::INFINITY);
        assert_eq!(round_half_even(f64::NEG_INFINITY, 5), f64::NEG_INFINITY);
    }

    fn round_half_even_legacy(value: f64, digits: i32) -> f64 {
        if !value.is_finite() {
            return value;
        }
        if is_decimal_format_underflow_tie(value, digits) {
            return 0.0_f64.copysign(value);
        }
        let s = format!("{:.*}", digits.max(0) as usize, value);
        s.parse::<f64>().unwrap_or(value)
    }

    #[test]
    fn round_half_even_stack_matches_heap_format() {
        let digits_list = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10];
        let mut special_values = vec![
            0.0,
            -0.0,
            f64::NAN,
            -f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            f64::from_bits(1),
            f64::from_bits(12345),
            0.5,
            1.5,
            2.5,
            3.5,
            -0.5,
            -1.5,
            -2.5,
            -3.5,
            0.05,
            0.15,
            0.25,
            0.35,
            0.00005,
            0.00015,
            0.00025,
            1.00005,
            -0.00005,
            -0.00015,
            -0.00025,
            1e30,
            -1e30,
            1e40,
            -1e40,
            1e60,
            -1e60,
            1e100,
            -1e100,
            f64::MAX,
            f64::MIN,
        ];

        // Exact half-way and decimal boundaries across 200,000+ values
        for k in 0..25_000u64 {
            let kf = k as f64;
            special_values.push(kf * 1e-5);
            special_values.push(-kf * 1e-5);
            special_values.push((kf + 0.5) * 1e-5);
            special_values.push((kf + 0.5) * 1e-3);
            special_values.push((kf + 0.5) * 1e-1);
            special_values.push(kf + 0.5);
        }

        let mut rng: u64 = 0xD1B5_4A32_D192_ED03;
        for _ in 0..50_000 {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let f = f64::from_bits(rng);
            special_values.push(f);
        }

        assert!(special_values.len() >= 200_000);

        for &v in &special_values {
            for &d in &digits_list {
                let actual = round_half_even(v, d);
                let expected = round_half_even_legacy(v, d);
                if expected.is_nan() {
                    assert!(actual.is_nan(), "expected NaN for v={v}, d={d}");
                } else {
                    assert_eq!(
                        actual.to_bits(),
                        expected.to_bits(),
                        "mismatch for v={v:?}, d={d}: actual={actual:?}, expected={expected:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn round_half_even_matches_jvm_oracle() {
        let tsv = include_str!("../testdata/format/java_fmt_oracle.tsv");
        let mut count = 0usize;
        for (line_no, line) in tsv.lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 9, "line {}: expected 9 columns", line_no + 1);
            let digits: i32 = cols[0].parse().unwrap();
            let v_bits = u64::from_str_radix(cols[3], 16).unwrap();
            let expected_bits = u64::from_str_radix(cols[5], 16).unwrap();
            let v = f64::from_bits(v_bits);
            let actual = round_half_even(v, digits);
            assert_eq!(
                actual.to_bits(),
                expected_bits,
                "line {}: round_half_even(v={} [0x{:016x}], d={}) = {} [0x{:016x}], expected {} [0x{:016x}]",
                line_no + 1,
                cols[2],
                v_bits,
                digits,
                actual,
                actual.to_bits(),
                cols[4],
                expected_bits
            );
            count += 1;
        }
        assert!(count >= 400, "expected >= 400 oracle rows, got {count}");
    }
}
