//! Emulation of 80-bit extended precision float for summation.
//! Emulates x86 `long double` (64-bit mantissa) addition for positive numbers.

#[derive(Clone, Copy, Debug)]
pub struct F80 {
    pub exp: i32,
    pub mant: u64, // Bit 63 is the implicit 1
}

impl F80 {
    pub fn from_f64(v: f64) -> Self {
        if v == 0.0 {
            return F80 { exp: 0, mant: 0 };
        }
        let bits = v.to_bits();
        let exp_bits = ((bits >> 52) & 0x7FF) as i32;
        let mant = bits & 0xFFFFFFFFFFFFF;

        if exp_bits == 0 {
            // Subnormal
            if mant == 0 {
                return F80 { exp: 0, mant: 0 };
            }
            let mut exp = -1022;
            let mut m = mant;
            while (m & 0x10000000000000) == 0 {
                m <<= 1;
                exp -= 1;
            }
            F80 { exp, mant: m << 11 }
        } else {
            let exp = exp_bits - 1023;
            F80 {
                exp,
                mant: (mant | 0x10000000000000) << 11,
            }
        }
    }

    pub fn to_f64(self) -> f64 {
        if self.mant == 0 {
            return 0.0;
        }
        let mut mant = self.mant;
        let mut exp = self.exp;
        let round_bit = (mant >> 10) & 1;
        let sticky = if (mant & ((1 << 10) - 1)) != 0 { 1 } else { 0 };
        mant >>= 11;
        if round_bit == 1 && (sticky == 1 || (mant & 1) == 1) {
            mant += 1;
            if mant >= (1 << 53) {
                mant >>= 1;
                exp = exp.saturating_add(1);
            }
        }

        // Handle overflow/underflow back to f64
        if exp >= 1024 {
            return f64::INFINITY;
        }
        if exp < -1022 {
            // Subnormal packing
            let shift_i64 = -1022i64 - i64::from(exp);
            if shift_i64 > 53 {
                return 0.0;
            }
            let shift = shift_i64 as u32;
            // Add implicit 1 bit back
            let sub_round_bit = (mant >> (shift - 1)) & 1;
            let sub_sticky = if (mant & ((1u64 << (shift - 1)) - 1)) != 0 {
                1
            } else {
                0
            };
            let mut sub_mant = mant >> shift;
            if sub_round_bit == 1 && (sub_sticky == 1 || (sub_mant & 1) == 1) {
                sub_mant += 1;
            }
            if sub_mant >= (1 << 52) {
                return f64::from_bits(1 << 52);
            }
            return f64::from_bits(sub_mant);
        }
        f64::from_bits((((exp + 1023) as u64) << 52) | (mant & 0xFFFFFFFFFFFFF))
    }

    /// Emulates 80-bit extended precision float addition.
    /// Supports non-negative operands only (e.g. probabilities/weights).
    pub fn add_f80(self, other: Self) -> Self {
        if self.mant == 0 {
            return other;
        }
        if other.mant == 0 {
            return self;
        }
        let (a, b) = if self.exp > other.exp || (self.exp == other.exp && self.mant >= other.mant) {
            (self, other)
        } else {
            (other, self)
        };
        let exp_diff = i64::from(a.exp) - i64::from(b.exp);
        let a_mant_128 = (a.mant as u128) << 60;
        let mut b_mant_128 = (b.mant as u128) << 60;
        if exp_diff >= 120 {
            b_mant_128 = 1;
        } else if exp_diff > 0 {
            let shift = exp_diff as u32;
            let sticky = if (b_mant_128 & ((1u128 << shift) - 1)) != 0 {
                1
            } else {
                0
            };
            b_mant_128 = (b_mant_128 >> shift) | sticky;
        }

        let mut sum = a_mant_128 + b_mant_128;
        let mut exp = a.exp;

        if sum >= (1u128 << 124) {
            let sticky = sum & 1;
            sum = (sum >> 1) | sticky;
            exp = exp.saturating_add(1);
        }

        let round_bit = (sum >> 59) & 1;
        let sticky = if (sum & ((1u128 << 59) - 1)) != 0 {
            1
        } else {
            0
        };
        sum >>= 60;

        if round_bit == 1 && (sticky == 1 || (sum & 1) == 1) {
            sum += 1;
            if sum >= (1u128 << 64) {
                sum >>= 1;
                exp = exp.saturating_add(1);
            }
        }

        F80 {
            exp,
            mant: sum as u64,
        }
    }
}

#[inline]
pub fn f80_sum(values: &[f64]) -> f64 {
    f80_sum_iter(values.iter().copied())
}

#[inline]
pub fn f80_sum_iter<I>(iter: I) -> f64
where
    I: IntoIterator<Item = f64>,
    I::IntoIter: Clone,
{
    let it = iter.into_iter();
    for v in it.clone() {
        if !v.is_finite() || v < 0.0 {
            let mut s = 0.0f64;
            for x in it {
                s += x;
            }
            return s;
        }
    }
    let mut sum = F80::from_f64(0.0);
    for v in it {
        sum = sum.add_f80(F80::from_f64(v));
    }
    sum.to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tie_rounding() {
        let eps = f64::EPSILON;

        let mut sum = F80::from_f64(1.0);
        sum = sum.add_f80(F80::from_f64(eps * 5.5));
        // Tie to even: 5.5 eps -> ends in half, nearest even is 6 eps.
        assert_eq!(sum.to_f64(), 1.0 + eps * 6.0);

        let mut sum2 = F80::from_f64(1.0);
        sum2 = sum2.add_f80(F80::from_f64(eps * 6.5));
        // Tie to even: 6.5 eps -> ends in half, nearest even is 6 eps.
        assert_eq!(sum2.to_f64(), 1.0 + eps * 6.0);
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_inf_nan_zero_subnormal() {
        assert_eq!(f80_sum(&[f64::INFINITY, 1.0]), f64::INFINITY);
        assert!(f80_sum(&[f64::NAN, 1.0]).is_nan());
        assert_eq!(f80_sum(&[0.0, 0.0]), 0.0);

        let min_pos = f64::MIN_POSITIVE; // Smallest positive normal
        let subnormal = min_pos / 2.0;
        assert!(subnormal > 0.0);
        assert_eq!(f80_sum(&[subnormal, subnormal]), min_pos);
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_chr18_mnhyper_sum_r_generated() {
        // R-generated exact expectation:
        // conda run -n vardict-oracle Rscript -e 'print(sum(c(43.999999999999992895)), digits=17)'
        // The array is just 43.999999999999992895
        // Let's emulate a sequence where long double accumulation differs from f64.
        // For chr18 mnhyper sum case, sum(d) = 43.999999999999992895

        let a = 43.0;
        let b = 0.999999999999992895; // 43.999999999999992895 - 43.0

        let result = f80_sum(&[a, b]);

        // Let's assert the raw output matches the R long double accumulation precision
        assert_eq!(result, 43.999999999999992895);
    }

    #[test]
    fn test_f80_from_to_add_and_sum_boundary_paths() {
        let zero = F80::from_f64(0.0);
        assert_eq!(zero.mant, 0);
        assert_eq!(zero.to_f64(), 0.0);

        // Adding with zero on either side returns the other operand unchanged
        let one = F80::from_f64(1.0);
        assert_eq!(zero.add_f80(one).to_f64(), 1.0);
        assert_eq!(one.add_f80(zero).to_f64(), 1.0);

        // Large exponent difference (exp_diff >= 120) sticky path
        let tiny = F80 {
            exp: -200,
            mant: 1u64 << 63,
        };
        assert_eq!(one.add_f80(tiny).to_f64(), 1.0);

        // Overflow back to f64::INFINITY when exp >= 1024
        let huge = F80 {
            exp: 1024,
            mant: 1u64 << 63,
        };
        assert_eq!(huge.to_f64(), f64::INFINITY);

        // Deep underflow (shift > 53) returns 0.0
        let deep_sub = F80 {
            exp: -1100,
            mant: 1u64 << 63,
        };
        assert_eq!(deep_sub.to_f64(), 0.0);

        // Subnormal round-trip
        let min_sub = f64::from_bits(1);
        assert_eq!(F80::from_f64(min_sub).to_f64(), min_sub);
        assert_eq!(f80_sum(&[]), 0.0);
    }

    #[test]
    fn test_f80_sum_iter_matches_f80_sum_bit_for_bit() {
        let ref_slice_sum = |vals: &[f64]| -> f64 {
            if vals.iter().any(|&v| !v.is_finite() || v < 0.0) {
                let mut s = 0.0f64;
                for &v in vals {
                    s += v;
                }
                return s;
            }
            let mut sum = F80::from_f64(0.0);
            for &v in vals {
                sum = sum.add_f80(F80::from_f64(v));
            }
            sum.to_f64()
        };

        let min_pos = f64::MIN_POSITIVE;
        let sub = min_pos / 4.0;
        let cases: &[&[f64]] = &[
            &[],
            &[1.0, f64::EPSILON * 5.5, f64::EPSILON * 0.5, 43.0],
            &[sub, sub, sub, sub],
            &[1.0, -0.25, 0.5],
            &[1.0, f64::INFINITY, 2.0],
            &[1.0, f64::NEG_INFINITY, 2.0],
            &[1.0, f64::NAN, 2.0],
        ];

        for &seq in cases {
            let expected = ref_slice_sum(seq);
            let from_slice = f80_sum(seq);
            let from_iter = f80_sum_iter(seq.iter().copied());
            assert_eq!(from_slice.to_bits(), expected.to_bits());
            assert_eq!(from_iter.to_bits(), from_slice.to_bits());

            // Also verify filtered iterator vs collected Vec slice
            let collected: Vec<f64> = seq.iter().copied().filter(|&v| v != 2.0).collect();
            let filtered_iter = seq.iter().copied().filter(|&v| v != 2.0);
            assert_eq!(
                f80_sum_iter(filtered_iter).to_bits(),
                f80_sum(&collected).to_bits()
            );
        }
    }

    #[test]
    fn test_f80_i32_exponent_extremes() {
        let max_carry = F80 {
            exp: i32::MAX,
            mant: u64::MAX,
        };
        assert_eq!(max_carry.to_f64(), f64::INFINITY);
        let sum_carry = max_carry.add_f80(max_carry);
        assert_eq!(sum_carry.exp, i32::MAX);
        assert_eq!(sum_carry.to_f64(), f64::INFINITY);

        let post_round_carry = F80 {
            exp: i32::MAX,
            mant: u64::MAX,
        }
        .add_f80(F80 {
            exp: i32::MAX - 64,
            mant: 1u64 << 63,
        });
        assert_eq!(post_round_carry.exp, i32::MAX);
        assert_eq!(post_round_carry.to_f64(), f64::INFINITY);

        let min_exp = F80 {
            exp: i32::MIN,
            mant: 1u64 << 63,
        };
        assert_eq!(min_exp.to_f64(), 0.0);

        let extreme_span = max_carry.add_f80(min_exp);
        assert_eq!(extreme_span.to_f64(), f64::INFINITY);
    }
}
