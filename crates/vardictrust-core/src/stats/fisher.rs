//! Fisher's Exact Test implementation mirroring VarDictJava's R-matching implementation.

use crate::types::round_half_even;

/// Computes Fisher's Exact Test p-values and conditional MLE odds ratios
/// exactly matching the Java reference codebase.
#[derive(Debug, Clone)]
pub struct FisherExact {
    m: i32,
    n: i32,
    k: i32,
    x: i32,
    lo: i32,
    hi: i32,
    logdc: Vec<f64>,
    p_value_less: f64,
    p_value_greater: f64,
    p_value_two_sided: f64,
    mu_at_1: f64,
    java_fisher: bool,
}

impl FisherExact {
    /// Constructs a new contingency table and performs calculations (R `fisher.test` mode).
    ///
    /// The contingency table has layout:
    /// ```text
    ///          Fwd       Rev
    /// Ref:   ref_fwd   ref_rev   |  m (total Ref)
    /// Alt:   alt_fwd   alt_rev   |  n (total Alt)
    ///        -----------------
    ///          k (Fwd)   (Rev)
    /// ```
    pub fn new(ref_fwd: i32, ref_rev: i32, alt_fwd: i32, alt_rev: i32) -> Self {
        Self::new_with_mode(ref_fwd, ref_rev, alt_fwd, alt_rev, false)
    }

    /// Constructs a new contingency table emulating `FisherExact.java` (`--fisher` mode),
    /// where `logdc` entries are rounded to 7 decimals (`roundHalfEven("0.0000000", value)`)
    /// and sums are accumulated in standard IEEE-754 `f64`.
    pub fn new_java(ref_fwd: i32, ref_rev: i32, alt_fwd: i32, alt_rev: i32) -> Self {
        Self::new_with_mode(ref_fwd, ref_rev, alt_fwd, alt_rev, true)
    }

    /// Constructs a new contingency table with explicit R (`false`) vs Java `--fisher` (`true`) mode.
    pub fn new_with_mode(
        ref_fwd: i32,
        ref_rev: i32,
        alt_fwd: i32,
        alt_rev: i32,
        java_fisher: bool,
    ) -> Self {
        let m = ref_fwd.saturating_add(ref_rev);
        let n = alt_fwd.saturating_add(alt_rev);
        let k = ref_fwd.saturating_add(alt_fwd);
        let x = ref_fwd;
        let lo = std::cmp::max(0, k.saturating_sub(n));
        let hi = std::cmp::min(k, m);

        let mut fisher = FisherExact {
            m,
            n,
            k,
            x,
            lo,
            hi,
            logdc: Vec::new(),
            p_value_less: 0.0,
            p_value_greater: 0.0,
            p_value_two_sided: 0.0,
            mu_at_1: f64::NAN,
            java_fisher,
        };

        fisher.logdc = fisher.logdc_dhyper(m, n, k);
        fisher.calculate_p_value();
        fisher
    }

    #[inline]
    fn support_len(&self) -> usize {
        if self.hi < self.lo {
            0
        } else {
            (self.hi as i64 - self.lo as i64 + 1) as usize
        }
    }

    /// Density of the central hypergeometric distribution on its support.
    fn logdc_dhyper(&self, m: i32, n: i32, k: i32) -> Vec<f64> {
        let len = self.support_len();
        let mut logdc = Vec::with_capacity(len);
        for i in 0..len {
            let element = self.lo.saturating_add(i as i32);
            if m + n == 0 {
                logdc.push(0.0);
                continue;
            }
            let mut value = if self.java_fisher {
                cm3_hypergeometric_log_probability(
                    element,
                    m.saturating_add(n),
                    m,
                    k,
                    self.lo,
                    self.hi,
                )
            } else {
                crate::stats::rmath::dhyper::dhyper(
                    element as f64,
                    m as f64,
                    n as f64,
                    k as f64,
                    true,
                )
            };
            if value.is_nan() {
                value = 0.0;
            }
            if self.java_fisher {
                value = round_half_even(value, 7);
            }
            logdc.push(value);
        }
        logdc
    }

    /// Determine the MLE for ncp by solving E(X) = x, where the expectation is with respect to H.
    fn mle(&self, x: f64) -> f64 {
        let eps = f64::EPSILON;
        if x == self.lo as f64 {
            return 0.0;
        }
        if x == self.hi as f64 {
            return f64::INFINITY;
        }
        // Function-local per-thread scratch, passed by `&mut` into every
        // `mnhyper` evaluation below (including all `uniroot` iterations).
        // Never shared across rayon workers; each call owns its buffer, so no
        // cross-thread aliasing is possible. MLE tolerance `eps.powf(0.25)`
        // is verbatim.
        let mut scratch = Vec::with_capacity(self.support_len());
        let mu = if x == self.x as f64 && !self.mu_at_1.is_nan() {
            self.mu_at_1
        } else {
            self.mnhyper_with(1.0, &mut scratch)
        };
        let tol = eps.powf(0.25);
        if mu > x {
            let f = |t: f64| self.mnhyper_with(t, &mut scratch) - x;
            crate::stats::zeroin::uniroot(0.0, 1.0, f, tol)
        } else if mu < x {
            let f = |t: f64| self.mnhyper_with(1.0 / t, &mut scratch) - x;
            1.0 / crate::stats::zeroin::uniroot(eps, 1.0, f, tol)
        } else {
            1.0
        }
    }

    /// Scratch-based mean evaluation. `scratch` is a function-local per-thread
    /// buffer (`&mut` only, never shared across threads); the products are
    /// formed in place in support order and summed with the same
    /// `java_stream_sum` vs `f80_sum` split as before, so results are
    /// element-for-element identical to the previous allocating version.
    pub(crate) fn mnhyper_with(&self, ncp: f64, scratch: &mut Vec<f64>) -> f64 {
        if ncp == 0.0 {
            return self.lo as f64;
        }
        if ncp.is_infinite() {
            return self.hi as f64;
        }
        self.dnhyper_into(ncp, scratch);
        for (i, prob) in scratch.iter_mut().enumerate() {
            let supp = self.lo.saturating_add(i as i32) as f64;
            *prob *= supp;
        }
        if self.java_fisher {
            java_stream_sum(scratch)
        } else {
            crate::stats::f80::f80_sum(scratch)
        }
    }

    /// Density of the non-central hypergeometric distribution on its support,
    /// written into the caller-owned scratch buffer `out` (function-local
    /// per-thread `&mut Vec`, never shared across rayon workers). The loop
    /// order, max-subtraction, summation choice, and normalization are
    /// identical to the previous allocating `dnhyper`; only the allocations
    /// are gone (single buffer reused in place).
    fn dnhyper_into(&self, ncp: f64, out: &mut Vec<f64>) {
        let len = self.support_len();
        out.clear();
        out.reserve(len);
        let ln_ncp = ncp.ln();
        for i in 0..len {
            let logdc_val = self.logdc.get(i).copied().unwrap_or(0.0);
            let supp_val = self.lo.saturating_add(i as i32) as f64;
            out.push(logdc_val + ln_ncp * supp_val);
        }
        let mut max_result = f64::NEG_INFINITY;
        for &val in out.iter() {
            if val > max_result {
                max_result = val;
            }
        }
        for val in out.iter_mut() {
            *val = (*val - max_result).exp();
        }
        let sum = if self.java_fisher {
            java_stream_sum(out)
        } else {
            crate::stats::f80::f80_sum(out)
        };
        for val in out.iter_mut() {
            *val /= sum;
        }
    }

    /// Returns the odds ratio formatted as String, matching Java (`--fisher`) or R outputs.
    /// Ports `FisherExact.java:149-156`.
    pub fn get_odds_ratio(&self) -> String {
        let odd_ratio = self.mle(self.x as f64);
        if odd_ratio.is_infinite() {
            "Inf".to_string()
        } else if odd_ratio == odd_ratio.round() {
            format!("{:.0}", odd_ratio)
        } else {
            let rounded = self.round_result(odd_ratio);
            if self.java_fisher {
                crate::variations::java_double_to_string(rounded)
            } else {
                format!("{}", rounded)
            }
        }
    }

    /// Returns the two-sided p-value.
    pub fn get_p_value(&self) -> f64 {
        self.round_result(self.p_value_two_sided)
    }

    /// Returns the one-sided "greater" p-value.
    pub fn get_p_value_greater(&self) -> f64 {
        self.round_result(self.p_value_greater)
    }

    /// Returns the one-sided "less" p-value.
    pub fn get_p_value_less(&self) -> f64 {
        self.round_result(self.p_value_less)
    }

    #[inline]
    fn round_result(&self, value: f64) -> f64 {
        if self.java_fisher {
            round_as_r_java(value)
        } else {
            round_as_r_fround5(value)
        }
    }

    fn calculate_p_value(&mut self) {
        self.p_value_less = self.pnhyper(self.x, false);
        self.p_value_greater = self.pnhyper(self.x, true);

        let rel_err = 1.0 + 1e-7;
        // Function-local per-thread scratch buffers, passed by `&mut` only
        // (never shared across rayon workers). The collect-then-sum shape is
        // kept verbatim: candidates are gathered in density order, then summed
        // with the same `iter().sum()` vs `f80_sum` split as before.
        let mut density = Vec::with_capacity(self.support_len());
        self.dnhyper_into(1.0, &mut density);
        let target_idx = self
            .x
            .checked_sub(self.lo)
            .and_then(|d| usize::try_from(d).ok())
            .unwrap_or(usize::MAX);
        let target = density.get(target_idx).copied().unwrap_or(0.0);
        let threshold = target * rel_err;
        if self.x == self.lo || self.x == self.hi {
            density.retain(|&el| el <= threshold);
            self.p_value_two_sided = if self.java_fisher {
                density.iter().copied().sum()
            } else {
                crate::stats::f80::f80_sum(&density)
            };
        } else {
            let mut sum_elements = Vec::with_capacity(density.len());
            for &el in &density {
                if el <= threshold {
                    sum_elements.push(el);
                }
            }
            self.p_value_two_sided = if self.java_fisher {
                sum_elements.iter().copied().sum()
            } else {
                crate::stats::f80::f80_sum(&sum_elements)
            };
            for (i, prob) in density.iter_mut().enumerate() {
                let supp = self.lo.saturating_add(i as i32) as f64;
                *prob *= supp;
            }
            self.mu_at_1 = if self.java_fisher {
                java_stream_sum(&density)
            } else {
                crate::stats::f80::f80_sum(&density)
            };
        }
    }

    fn pnhyper(&self, q: i32, upper_tail: bool) -> f64 {
        if self.m + self.n == 0 {
            return 1.0;
        }
        if self.java_fisher {
            let pop = self.m.saturating_add(self.n);
            let prob = |i: i32| -> f64 {
                let lp =
                    cm3_hypergeometric_log_probability(i, pop, self.m, self.k, self.lo, self.hi);
                if lp == f64::NEG_INFINITY {
                    0.0
                } else {
                    lp.exp()
                }
            };
            if upper_tail {
                if q <= self.lo {
                    return 1.0;
                }
                if q > self.hi {
                    return 0.0;
                }
                let mut x0 = self.hi;
                let mut ret = prob(x0);
                while x0 > q {
                    x0 -= 1;
                    ret += prob(x0);
                }
                return ret;
            } else {
                if q < self.lo {
                    return 0.0;
                }
                if q >= self.hi {
                    return 1.0;
                }
                let mut x0 = self.lo;
                let mut ret = prob(x0);
                while x0 < q {
                    x0 += 1;
                    ret += prob(x0);
                }
                return ret;
            }
        }
        if upper_tail {
            crate::stats::rmath::phyper::phyper(
                (q - 1) as f64,
                self.m as f64,
                self.n as f64,
                self.k as f64,
                false,
            )
        } else {
            crate::stats::rmath::phyper::phyper(
                q as f64,
                self.m as f64,
                self.n as f64,
                self.k as f64,
                true,
            )
        }
    }
}

/// Apache Commons Math 3.6.1 `SaddlePointExpansion.EXACT_STIRLING_ERRORS`.
#[allow(clippy::excessive_precision)]
const CM3_EXACT_STIRLING_ERRORS: [f64; 31] = [
    0.0,
    0.15342640972002736,
    0.08106146679532726,
    0.05481412105191765,
    0.0413406959554093,
    0.03316287351993629,
    0.02767792568499834,
    0.023746163656297496,
    0.020790672103765093,
    0.018488450532673187,
    0.016644691189821193,
    0.015134973221917378,
    0.013876128823070748,
    0.012810465242920227,
    0.01189670994589177,
    0.011104559758206917,
    0.010411265261972096,
    0.009799416126158804,
    0.009255462182712733,
    0.008768700134139386,
    0.00833056343336287,
    0.00793411456431402,
    0.007573675487951841,
    0.007244554301320383,
    0.00694284010720953,
    0.006665247032707682,
    0.006408994188004207,
    0.006171712263039458,
    0.0059513701127588475,
    0.0057462165130101155,
    0.005554733551962801,
];

/// Ports `org.apache.commons.math3.distribution.SaddlePointExpansion.getStirlingError`.
#[inline]
#[allow(clippy::excessive_precision)]
fn cm3_stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            return CM3_EXACT_STIRLING_ERRORS
                .get(z2 as usize)
                .copied()
                .unwrap_or(0.0);
        }
    }
    let z2 = z * z;
    (0.08333333333333333
        - (0.002777777777777778
            - (7.936507936507937e-4 - (5.952380952380953e-4 - 8.417508417508417e-4 / z2) / z2)
                / z2)
            / z2)
        / z
}

/// Ports `org.apache.commons.math3.distribution.SaddlePointExpansion.getDeviancePart`.
#[inline]
fn cm3_deviance_part(x: f64, mu: f64) -> f64 {
    if x == 0.0 && mu == 0.0 {
        return 0.0;
    }
    if (x - mu).abs() < 0.1 * (x + mu) {
        let d = x - mu;
        let v = d / (x + mu);
        let mut s1 = v * d;
        let mut ej = 2.0 * x * v;
        let v2 = v * v;
        for j in 1..1000i32 {
            let s = s1;
            ej *= v2;
            s1 = s + ej / (j * 2 + 1) as f64;
            if s1 == s {
                break;
            }
        }
        s1
    } else {
        x * (x / mu).ln() + mu - x
    }
}

/// Ports `org.apache.commons.math3.distribution.SaddlePointExpansion.logBinomialProbability`.
#[inline]
#[allow(clippy::excessive_precision)]
fn cm3_log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    let nf = n as f64;
    if x == 0 {
        if p < 0.1 {
            -cm3_deviance_part(nf, nf * q) - nf * p
        } else {
            nf * q.ln()
        }
    } else if x == n {
        if q < 0.1 {
            -cm3_deviance_part(nf, nf * p) - nf * q
        } else {
            nf * p.ln()
        }
    } else {
        let xf = x as f64;
        let n_minus_x = (n - x) as f64;
        let ret = cm3_stirling_error(nf)
            - cm3_stirling_error(xf)
            - cm3_stirling_error(n_minus_x)
            - cm3_deviance_part(xf, nf * p)
            - cm3_deviance_part(n_minus_x, nf * q);
        let f = std::f64::consts::TAU * xf * n_minus_x / nf;
        -0.5 * f.ln() + ret
    }
}

/// Ports `org.apache.commons.math3.distribution.HypergeometricDistribution.logProbability`.
#[inline]
fn cm3_hypergeometric_log_probability(
    x: i32,
    population_size: i32,
    number_of_successes: i32,
    sample_size: i32,
    lo: i32,
    hi: i32,
) -> f64 {
    if x < lo || x > hi {
        return f64::NEG_INFINITY;
    }
    let pop = population_size as f64;
    let p = (sample_size as f64) / pop;
    let q = ((population_size - sample_size) as f64) / pop;
    let p1 = cm3_log_binomial_probability(x, number_of_successes, p, q);
    let p2 =
        cm3_log_binomial_probability(sample_size - x, population_size - number_of_successes, p, q);
    let p3 = cm3_log_binomial_probability(sample_size, population_size, p, q);
    p1 + p2 - p3
}

/// Emulates OpenJDK's `DoubleStream.sum()` (`Collectors.sumWithCompensation` + `computeFinalSum`),
/// used by `FisherExact.java:112` and `FisherExact.java:129`.
fn java_stream_sum(values: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    let mut comp = 0.0_f64;
    let mut simple_sum = 0.0_f64;
    for &v in values {
        let tmp = v - comp;
        let velvel = sum + tmp;
        comp = (velvel - sum) - tmp;
        sum = velvel;
        simple_sum += v;
    }
    let tmp = sum - comp;
    if tmp.is_nan() && simple_sum.is_infinite() {
        simple_sum
    } else {
        tmp
    }
}

/// Ports `FisherExact.round_as_r` (`FisherExact.java:164-169`) for `--fisher` mode.
fn round_as_r_java(value: f64) -> f64 {
    let multiplied = value * 100000.0;
    let rounded = round_half_even(multiplied, 0) / 100000.0;
    if rounded == 0.0 {
        0.0
    } else if rounded == 1.0 {
        1.0
    } else {
        rounded
    }
}

/// Ports R 4.x `round(value, 5)` (`nmath/fround.c`) for non-`--fisher` R mode:
/// compares the two bounding 5-decimal candidates `floor(x * 1e5) / 1e5` and
/// `ceil(x * 1e5) / 1e5` directly against `x` in IEEE-754 `f64`, breaking exact
/// ties to the even integer significand.
fn round_as_r_fround5(value: f64) -> f64 {
    if !value.is_finite() || value.abs() >= (1u64 << 52) as f64 / 100000.0 {
        return round_as_r_java(value);
    }
    let scaled = value * 100000.0;
    let k_lo = scaled.floor();
    let k_hi = scaled.ceil();
    let x_lo = k_lo / 100000.0;
    let x_hi = k_hi / 100000.0;
    let d_lo = (value - x_lo).abs();
    let d_hi = (x_hi - value).abs();
    let rounded = if d_lo < d_hi {
        x_lo
    } else if d_hi < d_lo {
        x_hi
    } else if (k_lo as i64) & 1 == 0 {
        x_lo
    } else {
        x_hi
    };
    if rounded == 0.0 {
        0.0
    } else if rounded == 1.0 {
        1.0
    } else {
        rounded
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logdc_case1() {
        let fisher = FisherExact::new(11, 12, 1, 2);
        let logdc = &fisher.logdc;
        let expected = [-2.4696392, -1.0345547, -0.8675006, -1.9661129];
        assert_eq!(logdc.len(), expected.len());
        for (a, b) in logdc.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-6);
        }
    }

    #[test]
    fn test_s7_sbf_p_value() {
        let fisher = FisherExact::new(418, 421, 3, 2);
        assert_eq!(fisher.get_p_value(), 0.68601);
    }

    #[test]
    fn test_reference_cases() {
        let test_cases = vec![
            ((121, 55, 18, 23), 0.00378, 0.99908, 0.00287, 2.79657),
            ((121, 5, 18, 23), 0.0, 1.0, 0.0, 29.86184),
            ((37, 76, 1, 1), 1.0, 0.55362, 0.89275, 0.49015),
            ((0, 0, 1, 0), 1.0, 1.0, 1.0, 0.0),
            ((1, 0, 1, 0), 1.0, 1.0, 1.0, 0.0),
            ((0, 0, 0, 0), 1.0, 1.0, 1.0, 0.0),
            ((0, 1, 1, 0), 1.0, 0.5, 1.0, 0.0),
            ((1, 1, 1, 0), 1.0, 0.66667, 1.0, 0.0),
            ((1, 1, 1, 1), 1.0, 0.83333, 0.83333, 1.0),
            ((10, 10, 10, 1), 0.04722, 0.02599, 0.99802, 0.10703),
            ((10, 10, 10, 0), 0.01099, 0.00615, 1.0, 0.0),
            ((69, 1, 74, 95), 0.0, 1.0, 0.0, 87.68597),
            ((41, 86, 1, 1), 0.54688, 0.54687, 0.89571, 0.47973),
            ((130, 189, 1, 0), 0.40937, 0.40937, 1.0, 0.0),
            ((83, 40, 1, 2), 0.25746, 0.96473, 0.25746, 4.09908),
            ((74, 117, 1, 0), 0.39062, 0.39063, 1.0, 0.0),
            ((60, 62, 2, 0), 0.49593, 0.24797, 1.0, 0.0),
            ((43, 83, 1, 1), 1.0, 0.57111, 0.88361, 0.52091),
            ((78, 40, 1, 1), 1.0, 0.88515, 0.56849, 1.93844),
        ];

        for ((r_fwd, r_rev, a_fwd, a_rev), p_two, p_less, p_greater, odd_ratio) in test_cases {
            let fisher = FisherExact::new(r_fwd, r_rev, a_fwd, a_rev);
            let p_two_calc = fisher.get_p_value();
            let p_less_calc = fisher.get_p_value_less();
            let p_greater_calc = fisher.get_p_value_greater();
            let odd_ratio_calc: f64 = fisher.get_odds_ratio().parse().unwrap_or(0.0);

            // Compare rounded to 5 decimal places
            let r_p_two = round_half_even(p_two_calc, 5);
            let r_p_less = round_half_even(p_less_calc, 5);
            let r_p_greater = round_half_even(p_greater_calc, 5);
            let r_odd_ratio = round_half_even(odd_ratio_calc, 5);

            println!(
                "counts: ({}, {}, {}, {}) | raw_p_two: {}, raw_p_less: {}, raw_p_greater: {}, raw_odd: {}",
                r_fwd, r_rev, a_fwd, a_rev, fisher.p_value_two_sided, fisher.p_value_less, fisher.p_value_greater, fisher.mle(fisher.x as f64)
            );

            assert!(
                (r_p_two - p_two).abs() <= 1.1e-5,
                "Fail two-sided on ({}, {}, {}, {}): expected {}, got {}",
                r_fwd,
                r_rev,
                a_fwd,
                a_rev,
                p_two,
                r_p_two
            );
            assert!(
                (r_p_less - p_less).abs() <= 1.1e-5,
                "Fail less on ({}, {}, {}, {}): expected {}, got {}",
                r_fwd,
                r_rev,
                a_fwd,
                a_rev,
                p_less,
                r_p_less
            );
            assert!(
                (r_p_greater - p_greater).abs() <= 1.1e-5,
                "Fail greater on ({}, {}, {}, {}): expected {}, got {}",
                r_fwd,
                r_rev,
                a_fwd,
                a_rev,
                p_greater,
                r_p_greater
            );
            assert!(
                (r_odd_ratio - odd_ratio).abs() <= 1.1e-5,
                "Fail odd ratio on ({}, {}, {}, {}): expected {}, got {}",
                r_fwd,
                r_rev,
                a_fwd,
                a_rev,
                odd_ratio,
                r_odd_ratio
            );
        }
    }

    #[test]
    fn test_cross_check_fishers_exact() {
        use fishers_exact::fishers_exact;

        let mut count = 0;
        // Step by 4 to sample the space quickly and keep the test fast (<0.05s)
        for r_fwd in (0..35).step_by(4) {
            for r_rev in (0..35).step_by(4) {
                for a_fwd in (0..35).step_by(4) {
                    for a_rev in (0..35).step_by(4) {
                        if r_fwd + r_rev + a_fwd + a_rev == 0 {
                            continue;
                        }
                        let fisher = FisherExact::new(r_fwd, r_rev, a_fwd, a_rev);
                        if let Ok(p) =
                            fishers_exact(&[r_fwd as u32, r_rev as u32, a_fwd as u32, a_rev as u32])
                        {
                            let our_p = fisher.p_value_two_sided;
                            let our_p_less = fisher.p_value_less;
                            let our_p_greater = fisher.p_value_greater;

                            assert!(
                                (our_p - p.two_tail_pvalue).abs() <= 1e-4,
                                "Mismatch two-tail at ({}, {}, {}, {}): ours={}, crate={}",
                                r_fwd,
                                r_rev,
                                a_fwd,
                                a_rev,
                                our_p,
                                p.two_tail_pvalue
                            );
                            assert!(
                                (our_p_less - p.less_pvalue).abs() <= 1e-4,
                                "Mismatch less at ({}, {}, {}, {}): ours={}, crate={}",
                                r_fwd,
                                r_rev,
                                a_fwd,
                                a_rev,
                                our_p_less,
                                p.less_pvalue
                            );
                            assert!(
                                (our_p_greater - p.greater_pvalue).abs() <= 1e-4,
                                "Mismatch greater at ({}, {}, {}, {}): ours={}, crate={}",
                                r_fwd,
                                r_rev,
                                a_fwd,
                                a_rev,
                                our_p_greater,
                                p.greater_pvalue
                            );
                            count += 1;
                        }
                    }
                }
            }
        }
        println!("Cross-checked {} cases with fishers_exact crate", count);
    }

    #[test]
    fn test_chr10_odds_ratio() {
        let f = FisherExact::new(27, 17, 3, 8);
        assert_eq!(f.get_odds_ratio(), "4.12268");
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_chrx_sbf() {
        let f = FisherExact::new(86, 41, 0, 1);
        let p = f.p_value_two_sided;

        // The exact raw value generated by R (and now correctly reproduced by Rust via dhyper port)
        let r_raw = 0.32812500000000006;

        assert!(
            (p - r_raw).abs() < 1e-15,
            "Raw p-value {} deviated too far from exact R raw {}",
            p,
            r_raw
        );

        // The formatted value must be exactly 0.32813 (banker's rounding of 0.32812500000000006)
        assert_eq!(f.get_p_value(), 0.32813);
    }

    #[test]
    fn test_fisher_chr1_66987255() {
        let f = FisherExact::new(0, 42, 2, 85);
        assert_eq!(f.get_p_value_less(), 0.45313);
    }

    #[test]
    fn test_fisher_chr8_132972774() {
        let f = FisherExact::new(7, 8, 1, 4);
        assert_eq!(f.get_odds_ratio(), "3.30211");
    }

    #[test]
    fn test_fisher_chr16_16364166_sor() {
        let f = FisherExact::new(5, 21, 1, 12);
        assert_eq!(f.get_odds_ratio(), "2.79182");
    }

    #[test]
    fn test_fisher_chr17_76093187_sor() {
        let f = FisherExact::new(1, 65, 2, 31);
        assert_eq!(f.get_odds_ratio(), "0.24227");
    }

    #[test]
    fn test_java_fisher_ultradeep_parity() {
        let odd_cases = [
            ((3506, 936, 2, 2), "3.74441"),
            ((3463, 353, 1, 1), "9.79777"),
            ((1327, 50, 1, 1), "26.2886"),
            ((8742, 9897, 2, 6), "2.64965"),
            ((3150, 4037, 1, 1), "0.78032"),
            ((7143, 8187, 1, 2), "1.74489"),
            ((6051, 7004, 1, 2), "1.72779"),
            ((11019, 7350, 2, 12), "8.99336"),
            ((6489, 5960, 1, 5), "5.44288"),
            ((4617, 4347, 1, 4), "4.24747"),
            ((4636, 4638, 2, 16), "7.99514"),
            ((4550, 3886, 1, 7), "8.19429"),
            ((3742, 2601, 1, 4), "5.75367"),
            ((3601, 2577, 1, 4), "5.58853"),
            ((2530, 3095, 1, 4), "3.26911"),
            ((3931, 6596, 1, 4), "2.38369"),
            ((1874, 194, 1, 1), "9.63757"),
            ((777, 572, 1, 4), "5.42793"),
            ((1792, 2602, 1, 3), "2.06578"),
            ((10, 8, 1, 3), "3.53775"),
            ((11, 1, 9, 3), "3.4803"),
            ((11, 2, 7, 6), "4.43009"),
            ((18, 14, 2, 6), "3.73129"),
        ];
        for ((rf, rr, af, ar), expected_or) in odd_cases {
            let f = FisherExact::new_java(rf, rr, af, ar);
            assert_eq!(
                f.get_odds_ratio(),
                expected_or,
                "failed odds ratio for ({rf}, {rr}, {af}, {ar})"
            );
        }

        let pval_cases = [((2437, 3314, 0, 3), 0.26742), ((730, 497, 1, 2), 0.56967)];
        for ((rf, rr, af, ar), expected_p) in pval_cases {
            let f = FisherExact::new_java(rf, rr, af, ar);
            assert_eq!(
                f.get_p_value(),
                expected_p,
                "failed p-value for ({rf}, {rr}, {af}, {ar})"
            );
        }
    }

    #[test]
    fn linux_x86_64_r_fisher_canaries() {
        assert_eq!(FisherExact::new(9, 6, 5, 1).get_odds_ratio(), "0.31603");
        assert_eq!(FisherExact::new(27, 17, 3, 8).get_odds_ratio(), "4.12268");
        assert_eq!(FisherExact::new(10, 8, 1, 3).get_odds_ratio(), "3.53774");
        assert_eq!(FisherExact::new(2, 150, 7, 183).get_odds_ratio(), "0.34953");
    }

    #[test]
    fn fisher_java_odds_ratio_double_to_string_pins() {
        assert_eq!(
            FisherExact::new_java(2, 300, 300, 2).get_odds_ratio(),
            "1.0E-4"
        );
        assert_eq!(
            FisherExact::new_java(2, 400, 400, 2).get_odds_ratio(),
            "6.0E-5"
        );
        assert_eq!(
            FisherExact::new_java(1, 2, 1200, 1).get_odds_ratio(),
            "6.2E-4"
        );
        assert_eq!(
            FisherExact::new_java(1, 100, 100, 1).get_odds_ratio(),
            "1.6E-4"
        );
        assert_eq!(FisherExact::new_java(1, 1, 1, 8).get_odds_ratio(), "6.0");
        assert_eq!(FisherExact::new_java(8, 1, 1, 1).get_odds_ratio(), "6.0");
        assert_eq!(
            FisherExact::new_java(200, 1, 4, 800).get_odds_ratio(),
            "16384.0"
        );
        assert_eq!(
            FisherExact::new_java(300, 3, 3, 1200).get_odds_ratio(),
            "16384.0"
        );
        // R mode keeps non-scientific formatting ("0.0001")
        assert_eq!(FisherExact::new(2, 300, 300, 2).get_odds_ratio(), "0.0001");
    }

    #[test]
    fn fisher_java_oracle_table() {
        let tsv = include_str!("../../testdata/fisher/java_fisher_oracle.tsv");
        let mut count = 0usize;
        for (line_no, line) in tsv.lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 11, "line {}: expected 11 columns", line_no + 1);
            let rf: i32 = cols[0].parse().unwrap();
            let rr: i32 = cols[1].parse().unwrap();
            let af: i32 = cols[2].parse().unwrap();
            let ar: i32 = cols[3].parse().unwrap();
            let p_two_bits = u64::from_str_radix(cols[5], 16).unwrap();
            let p_less_bits = u64::from_str_radix(cols[7], 16).unwrap();
            let p_greater_bits = u64::from_str_radix(cols[9], 16).unwrap();
            let expected_or = cols[10];

            let fe = FisherExact::new_java(rf, rr, af, ar);
            assert_eq!(
                fe.get_p_value().to_bits(),
                p_two_bits,
                "line {} ({rf},{rr},{af},{ar}) p_two: got {} [0x{:016x}], expected {} [0x{:016x}]",
                line_no + 1,
                fe.get_p_value(),
                fe.get_p_value().to_bits(),
                cols[4],
                p_two_bits
            );
            assert_eq!(
                fe.get_p_value_less().to_bits(),
                p_less_bits,
                "line {} ({rf},{rr},{af},{ar}) p_less: got {} [0x{:016x}], expected {} [0x{:016x}]",
                line_no + 1,
                fe.get_p_value_less(),
                fe.get_p_value_less().to_bits(),
                cols[6],
                p_less_bits
            );
            assert_eq!(
                fe.get_p_value_greater().to_bits(),
                p_greater_bits,
                "line {} ({rf},{rr},{af},{ar}) p_greater: got {} [0x{:016x}], expected {} [0x{:016x}]",
                line_no + 1,
                fe.get_p_value_greater(),
                fe.get_p_value_greater().to_bits(),
                cols[8],
                p_greater_bits
            );
            assert_eq!(
                fe.get_odds_ratio(),
                expected_or,
                "line {} ({rf},{rr},{af},{ar}) odd_ratio: got {}, expected {}",
                line_no + 1,
                fe.get_odds_ratio(),
                expected_or
            );
            count += 1;
        }
        assert!(count >= 250, "expected >= 250 oracle rows, got {count}");
    }

    #[test]
    fn fisher_r_oracle_table() {
        let tsv = include_str!("../../testdata/fisher/r_fisher_oracle.tsv");
        let mut count = 0usize;
        for (line_no, line) in tsv.lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 8, "line {}: expected 8 columns", line_no + 1);
            let rf: i32 = cols[0].parse().unwrap();
            let rr: i32 = cols[1].parse().unwrap();
            let af: i32 = cols[2].parse().unwrap();
            let ar: i32 = cols[3].parse().unwrap();
            let expected_p_two: f64 = cols[4].parse().unwrap();
            let expected_p_less: f64 = cols[5].parse().unwrap();
            let expected_p_greater: f64 = cols[6].parse().unwrap();
            let expected_or = cols[7];

            let fe = FisherExact::new(rf, rr, af, ar);
            assert_eq!(
                fe.get_p_value().to_bits(),
                expected_p_two.to_bits(),
                "line {} ({rf},{rr},{af},{ar}) R p_two: got {}, expected {}",
                line_no + 1,
                fe.get_p_value(),
                expected_p_two
            );
            assert_eq!(
                fe.get_p_value_less().to_bits(),
                expected_p_less.to_bits(),
                "line {} ({rf},{rr},{af},{ar}) R p_less: got {}, expected {}",
                line_no + 1,
                fe.get_p_value_less(),
                expected_p_less
            );
            assert_eq!(
                fe.get_p_value_greater().to_bits(),
                expected_p_greater.to_bits(),
                "line {} ({rf},{rr},{af},{ar}) R p_greater: got {}, expected {}",
                line_no + 1,
                fe.get_p_value_greater(),
                expected_p_greater
            );
            assert_eq!(
                fe.get_odds_ratio(),
                expected_or,
                "line {} ({rf},{rr},{af},{ar}) R odd_ratio: got {}, expected {}",
                line_no + 1,
                fe.get_odds_ratio(),
                expected_or
            );
            count += 1;
        }
        assert!(count >= 250, "expected >= 250 oracle rows, got {count}");
    }

    #[test]
    fn mu_at_1_reuse_and_boundary_invariants() {
        let interior_cases = [
            (
                (121, 55, 18, 23),
                (0.00378, 0.99908, 0.00287, "2.79657"),
                (0.00378, 0.99908, 0.00287, "2.79657"),
            ),
            (
                (83, 40, 1, 2),
                (0.25746, 0.96473, 0.25746, "4.09908"),
                (0.25746, 0.96473, 0.25746, "4.09908"),
            ),
            (
                (10, 10, 10, 1),
                (0.04722, 0.02599, 0.99802, "0.10703"),
                (0.04722, 0.02599, 0.99802, "0.10703"),
            ),
        ];

        let mut scratch = Vec::new();
        for ((rf, rr, af, ar), expected_r, expected_java) in interior_cases {
            for (java_mode, (exp_p, exp_less, exp_greater, exp_or)) in
                [(false, expected_r), (true, expected_java)]
            {
                let f = if java_mode {
                    FisherExact::new_java(rf, rr, af, ar)
                } else {
                    FisherExact::new(rf, rr, af, ar)
                };
                assert!(f.x > f.lo && f.x < f.hi);
                assert!(!f.mu_at_1.is_nan());
                assert_eq!(
                    f.mu_at_1.to_bits(),
                    f.mnhyper_with(1.0, &mut scratch).to_bits()
                );
                assert_eq!(f.get_p_value(), exp_p);
                assert_eq!(f.get_p_value_less(), exp_less);
                assert_eq!(f.get_p_value_greater(), exp_greater);
                assert_eq!(f.get_odds_ratio(), exp_or);

                let custom_x = if f.x - 1 > f.lo {
                    (f.x - 1) as f64
                } else {
                    (f.x + 1) as f64
                };
                assert_ne!(custom_x, f.x as f64);
                let mut f_uncached = f.clone();
                f_uncached.mu_at_1 = f64::NAN;
                assert_eq!(
                    f.mle(f.x as f64).to_bits(),
                    f_uncached.mle(f.x as f64).to_bits()
                );
                let mut f_poisoned = f.clone();
                f_poisoned.mu_at_1 = 1.0e300;
                assert_eq!(
                    f.mle(custom_x).to_bits(),
                    f_poisoned.mle(custom_x).to_bits()
                );
            }
        }

        let boundary_cases = [
            ((0, 0, 0, 0), (1.0, 1.0, 1.0, "0"), (1.0, 1.0, 1.0, "0")),
            (
                (0, 10, 10, 10),
                (0.01099, 0.00615, 1.0, "0"),
                (0.01099, 0.00615, 1.0, "0"),
            ),
            (
                (10, 0, 10, 10),
                (0.01099, 1.0, 0.00615, "Inf"),
                (0.01099, 1.0, 0.00615, "Inf"),
            ),
            (
                (10, 10, 10, 0),
                (0.01099, 0.00615, 1.0, "0"),
                (0.01099, 0.00615, 1.0, "0"),
            ),
        ];

        for ((rf, rr, af, ar), expected_r, expected_java) in boundary_cases {
            for (java_mode, (exp_p, exp_less, exp_greater, exp_or)) in
                [(false, expected_r), (true, expected_java)]
            {
                let f = if java_mode {
                    FisherExact::new_java(rf, rr, af, ar)
                } else {
                    FisherExact::new(rf, rr, af, ar)
                };
                assert!(f.x == f.lo || f.x == f.hi);
                assert!(f.mu_at_1.is_nan());
                assert_eq!(f.get_p_value(), exp_p);
                assert_eq!(f.get_p_value_less(), exp_less);
                assert_eq!(f.get_p_value_greater(), exp_greater);
                assert_eq!(f.get_odds_ratio(), exp_or);

                if f.hi - f.lo >= 2 {
                    let custom_x = (f.lo + 1) as f64;
                    assert_ne!(custom_x, f.x as f64);
                    let f_interior = if java_mode {
                        FisherExact::new_java(
                            f.lo + 1,
                            f.m - (f.lo + 1),
                            f.k - (f.lo + 1),
                            f.n - (f.k - (f.lo + 1)),
                        )
                    } else {
                        FisherExact::new(
                            f.lo + 1,
                            f.m - (f.lo + 1),
                            f.k - (f.lo + 1),
                            f.n - (f.k - (f.lo + 1)),
                        )
                    };
                    assert_eq!(
                        f.mle(custom_x).to_bits(),
                        f_interior.mle(custom_x).to_bits()
                    );
                }
            }
        }
    }
}
