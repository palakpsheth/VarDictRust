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

/// Banker's rounding to `digits` decimal places.
///
/// Matches Java's `DecimalFormat` with `RoundingMode.HALF_EVEN`.
/// This is critical for float-field equivalence with VarDictJava.
pub fn round_half_even(value: f64, digits: i32) -> f64 {
    if !value.is_finite() {
        return value;
    }
    let s = format!("{:.*}", digits.max(0) as usize, value);
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
}
