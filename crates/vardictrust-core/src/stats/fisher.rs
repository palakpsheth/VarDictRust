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
    support: Vec<i32>,
    logdc: Vec<f64>,
    p_value_less: f64,
    p_value_greater: f64,
    p_value_two_sided: f64,
}

impl FisherExact {
    /// Constructs a new contingency table and performs calculations.
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
        let m = ref_fwd + ref_rev;
        let n = alt_fwd + alt_rev;
        let k = ref_fwd + alt_fwd;
        let x = ref_fwd;
        let lo = std::cmp::max(0, k - n);
        let hi = std::cmp::min(k, m);
        let mut support = Vec::new();
        for j in lo..=hi {
            support.push(j);
        }

        let mut fisher = FisherExact {
            m,
            n,
            k,
            x,
            lo,
            hi,
            support,
            logdc: Vec::new(),
            p_value_less: 0.0,
            p_value_greater: 0.0,
            p_value_two_sided: 0.0,
        };

        fisher.logdc = fisher.logdc_dhyper(m, n, k);
        fisher.calculate_p_value();
        fisher
    }

    /// Density of the central hypergeometric distribution on its support.
    fn logdc_dhyper(&self, m: i32, n: i32, k: i32) -> Vec<f64> {
        let mut logdc = Vec::with_capacity(self.support.len());
        for &element in &self.support {
            if m + n == 0 {
                logdc.push(0.0);
                continue;
            }
            let mut value = crate::stats::rmath::dhyper::dhyper(
                element as f64,
                m as f64,
                n as f64,
                k as f64,
                true,
            );
            if value.is_nan() {
                value = 0.0;
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
        let mu = self.mnhyper(1.0);
        let root;
        let tol = eps.powf(0.25);
        if mu > x {
            let f = |t: f64| self.mnhyper(t) - x;
            root = crate::stats::zeroin::uniroot(0.0, 1.0, f, tol);
        } else if mu < x {
            let f = |t: f64| self.mnhyper(1.0 / t) - x;
            root = 1.0 / crate::stats::zeroin::uniroot(eps, 1.0, f, tol);
        } else {
            root = 1.0;
        }
        root
    }

    pub(crate) fn mnhyper(&self, ncp: f64) -> f64 {
        if ncp == 0.0 {
            return self.lo as f64;
        }
        if ncp.is_infinite() {
            return self.hi as f64;
        }
        let dnhyper_result = self.dnhyper(ncp);
        let mut supp_d = Vec::with_capacity(self.support.len());
        for (i, &supp) in self.support.iter().enumerate() {
            supp_d.push((supp as f64) * dnhyper_result.get(i).copied().unwrap_or(0.0));
        }
        crate::stats::f80::f80_sum(&supp_d)
    }

    fn dnhyper(&self, ncp: f64) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.support.len());
        let ln_ncp = ncp.ln();
        for i in 0..self.support.len() {
            let logdc_val = self.logdc.get(i).copied().unwrap_or(0.0);
            let supp_val = self.support.get(i).copied().unwrap_or(0) as f64;
            result.push(logdc_val + ln_ncp * supp_val);
        }
        let mut max_result = f64::NEG_INFINITY;
        for &val in &result {
            if val > max_result {
                max_result = val;
            }
        }
        let mut exponent_result = Vec::with_capacity(result.len());
        for val in result {
            let exp_val = (val - max_result).exp();
            exponent_result.push(exp_val);
        }
        let sum = crate::stats::f80::f80_sum(&exponent_result);
        let mut final_result = Vec::with_capacity(exponent_result.len());
        for val in exponent_result {
            final_result.push(val / sum);
        }
        final_result
    }

    /// Returns the odds ratio formatted as String, matching R outputs.
    pub fn get_odds_ratio(&self) -> String {
        let odd_ratio = self.mle(self.x as f64);
        if odd_ratio.is_infinite() {
            "Inf".to_string()
        } else if odd_ratio == odd_ratio.round() {
            format!("{:.0}", odd_ratio)
        } else {
            format!("{}", round_as_r(odd_ratio))
        }
    }

    /// Returns the two-sided p-value.
    pub fn get_p_value(&self) -> f64 {
        round_as_r(self.p_value_two_sided)
    }

    /// Returns the one-sided "greater" p-value.
    pub fn get_p_value_greater(&self) -> f64 {
        round_as_r(self.p_value_greater)
    }

    /// Returns the one-sided "less" p-value.
    pub fn get_p_value_less(&self) -> f64 {
        round_as_r(self.p_value_less)
    }

    fn calculate_p_value(&mut self) {
        self.p_value_less = self.pnhyper(self.x, false);
        self.p_value_greater = self.pnhyper(self.x, true);

        let rel_err = 1.0 + 1e-7;
        let d = self.dnhyper(1.0);
        let target_idx = (self.x - self.lo) as usize;
        let target = d.get(target_idx).copied().unwrap_or(0.0);
        let threshold = target * rel_err;
        let mut sum_elements = Vec::new();
        for el in d {
            if el <= threshold {
                sum_elements.push(el);
            }
        }
        self.p_value_two_sided = crate::stats::f80::f80_sum(&sum_elements);
    }

    fn pnhyper(&self, q: i32, upper_tail: bool) -> f64 {
        if self.m + self.n == 0 {
            return 1.0;
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

fn round_as_r(value: f64) -> f64 {
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
}
