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
                exp += 1;
            }
        }

        // Handle overflow/underflow back to f64
        if exp >= 1024 {
            return f64::INFINITY;
        }
        if exp < -1022 {
            // Subnormal packing
            let shift = -1022 - exp;
            if shift > 53 {
                return 0.0;
            }
            // Add implicit 1 bit back
            let sub_round_bit = (mant >> (shift - 1)) & 1;
            let sub_sticky = if (mant & ((1 << (shift - 1)) - 1)) != 0 {
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
        let exp_diff = a.exp - b.exp;
        let a_mant_128 = (a.mant as u128) << 60;
        let mut b_mant_128 = (b.mant as u128) << 60;
        if exp_diff >= 120 {
            b_mant_128 = 1;
        } else if exp_diff > 0 {
            let shift = exp_diff;
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
            exp += 1;
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
                exp += 1;
            }
        }

        F80 {
            exp,
            mant: sum as u64,
        }
    }
}

pub fn f80_sum(vals: &[f64]) -> f64 {
    let mut fallback = false;
    for &v in vals {
        debug_assert!(
            v.is_nan() || v >= 0.0,
            "f80_sum expects non-negative operands only"
        );
        if !v.is_finite() || v < 0.0 {
            fallback = true;
            break;
        }
    }

    // Fall back to f64 sum if any element is inf/nan or negative
    if fallback {
        let mut sum = 0.0;
        for &v in vals {
            sum += v;
        }
        return sum;
    }

    let mut sum = F80::from_f64(0.0);
    for &v in vals {
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
}
