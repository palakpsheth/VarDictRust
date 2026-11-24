//! Cleanroom `MIT` saddlepoint binomial and hypergeometric probability engine.
//!
//! Derived directly from published mathematical identities:
//! - Catherine Loader (2000), *"Fast and Accurate Computation of Binomial Probabilities"*,
//!   HP Laboratories Technical Report (`Eqs. 2, 5, 6`, and `§5` three-binomial
//!   hypergeometric representation).
//! - Euler-Maclaurin asymptotic expansion for Stirling's remainder
//!   $\delta(z) = \ln \Gamma(z + 1) - (z + \tfrac{1}{2})\ln z + z - \tfrac{1}{2}\ln(2\pi)$
//!   using Bernoulli numbers $B_{2m}$ ($c_m = \frac{B_{2m}}{2m(2m - 1)}$) and the exact
//!   Gamma functional equation $\Gamma(z + 1) = z \Gamma(z)$.
//! - Milton Abramowitz & Irene A. Stegun (1964), *Handbook of Mathematical Functions*,
//!   `§26.3.29` hypergeometric probability ratio recurrence.

/// Exact IEEE-754 `f64` nearest value of $\ln(2\pi) = 1.83787706640934548356\dots$
/// (avoiding the 1-ULP double-rounding loss of `std::f64::consts::TAU.ln()`).
const LN_2PI: f64 = 1.8378770664093456;

/// Exact IEEE-754 `f64` values of $\delta(k / 2)$ for $k = 0..=30$ ($n \in \{0, 0.5, 1.0, \dots, 15.0\}$),
/// where $\delta(z) = \ln \Gamma(z + 1) - (z + \tfrac{1}{2})\ln z + z - \tfrac{1}{2}\ln(2\pi)$ and $\delta(0) = 0$.
const STIRLING_HALF_INTEGER_ERRORS: [f64; 31] = [
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

// Euler-Maclaurin Bernoulli coefficients c_m = B_{2m} / (2m * (2m - 1)) for m = 1..=8.
const C1: f64 = 1.0 / 12.0;
const C2: f64 = -1.0 / 360.0;
const C3: f64 = 1.0 / 1260.0;
const C4: f64 = -1.0 / 1680.0;
const C5: f64 = 1.0 / 1188.0;
const C6: f64 = -691.0 / 360_360.0;
const C7: f64 = 1.0 / 156.0;
const C8: f64 = -3617.0 / 122_400.0;

#[inline]
fn prob_zero(log_p: bool) -> f64 {
    if log_p {
        f64::NEG_INFINITY
    } else {
        0.0
    }
}

#[inline]
fn prob_one(log_p: bool) -> f64 {
    if log_p {
        0.0
    } else {
        1.0
    }
}

/// Evaluates Stirling's log-factorial error
/// $\delta(n) = \ln \Gamma(n + 1) - (n + \tfrac{1}{2})\ln n + n - \tfrac{1}{2}\ln(2\pi)$
/// with $\delta(0) = 0$.
#[inline]
pub fn stirling_error(n: f64) -> f64 {
    if n.is_nan() || n < 0.0 {
        return f64::NAN;
    }
    if !n.is_finite() {
        return 0.0;
    }
    let twice = n * 2.0;
    if n <= 15.0 && twice == twice.trunc() {
        return STIRLING_HALF_INTEGER_ERRORS
            .get(twice as usize)
            .copied()
            .unwrap_or(0.0);
    }
    if n > 15.0 {
        if n > 2.0e7 {
            return C1 / n;
        }
        let inv_n2 = 1.0 / (n * n);
        let poly = C1
            + inv_n2
                * (C2
                    + inv_n2
                        * (C3
                            + inv_n2
                                * (C4
                                    + inv_n2
                                        * (C5 + inv_n2 * (C6 + inv_n2 * (C7 + inv_n2 * C8))))));
        return poly / n;
    }
    // Non-half-integer 0.0 < n < 15.0: exact Gamma(z + 1) = z * Gamma(z) shift by 16:
    // Gamma(n + 17) = Gamma(n + 1) * prod_{k=1}^{16} (n + k).
    let m = n + 16.0;
    let mut log_prod = 0.0;
    for k in 1..=16_u8 {
        log_prod += (n + f64::from(k)).ln();
    }
    stirling_error(m) + (m + 0.5) * m.ln() - (n + 0.5) * n.ln() - 16.0 - log_prod
}

/// Evaluates Loader's (2000, Eq. 5) deviance term $D_0(x, \mu) = x \ln(x / \mu) + \mu - x$.
#[inline]
pub fn binomial_deviance(x: f64, mu: f64) -> f64 {
    if !x.is_finite() || !mu.is_finite() || mu <= 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return mu;
    }
    if x == mu {
        return 0.0;
    }
    let sum_x_mu = x + mu;
    let near_diag = if sum_x_mu.is_infinite() {
        (0.25 * x - 0.25 * mu).abs() < 0.1 * (0.25 * x + 0.25 * mu)
    } else {
        (x - mu).abs() < 0.1 * sum_x_mu
    };
    if near_diag {
        let d = x - mu;
        let v = if sum_x_mu.is_infinite() {
            (0.25 * x - 0.25 * mu) / (0.25 * x + 0.25 * mu)
        } else {
            d / sum_x_mu
        };
        let mut s = d * v;
        if s == 0.0 {
            return 0.0;
        }
        let mut ej = if sum_x_mu.is_infinite() {
            (x * v) * 2.0
        } else {
            2.0 * x * v
        };
        let v2 = v * v;
        for j in 1..1000_u32 {
            ej *= v2;
            let s_next = s + ej / f64::from(2 * j + 1);
            if s_next == s {
                return s_next;
            }
            s = s_next;
        }
        return s;
    }
    let ratio = x / mu;
    let lg = if ratio.is_finite() && x > 0.0 {
        ratio.ln()
    } else {
        x.ln() - mu.ln()
    };
    if x > mu {
        x * (lg - 1.0) + mu
    } else {
        x * lg + mu - x
    }
}

/// Evaluates $(1 + x)^y$ accurately for $|x| \ll 1$ and small integer exponents $y \in \{0, 1, 2, 3, 4\}$.
///
/// For $|x| \le \tfrac{1}{2}$, decomposes $1 + x = u + \delta$ exactly via Dekker's `FastTwoSum`
/// ($u = \operatorname{fl}(1 + x)$, $\delta = x - (u - 1)$) so $(1 + x)^y = u^y \exp(y \ln(1 + \delta / u))$
/// preserves both exact binary bases and sub-ULP residuals when $|x| \ll 1$.
#[inline]
pub fn one_plus_pow(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if y == 0.0 || (x == 0.0 && y.is_finite()) {
        return 1.0;
    }
    if x < -1.0 {
        return (1.0 + x).powf(y);
    }
    if x == -1.0 {
        return 0.0_f64.powf(y);
    }
    if y == 1.0 {
        return 1.0 + x;
    }
    if y == 2.0 {
        return if x.abs() <= 0.5 {
            1.0 + x * (2.0 + x)
        } else {
            let u = 1.0 + x;
            u * u
        };
    }
    if y == 3.0 {
        return if x.abs() <= 0.5 {
            1.0 + x * (3.0 + x * (3.0 + x))
        } else {
            let u = 1.0 + x;
            u * u * u
        };
    }
    if y == 4.0 {
        return if x.abs() <= 0.5 {
            1.0 + x * (4.0 + x * (6.0 + x * (4.0 + x)))
        } else {
            let u = 1.0 + x;
            let u2 = u * u;
            u2 * u2
        };
    }
    if x.abs() <= 0.5 {
        let u = 1.0 + x;
        let delta = x - (u - 1.0);
        let base = u.powf(y);
        if delta == 0.0 {
            return base;
        }
        if base == 0.0 || !base.is_finite() {
            return (y * x.ln_1p()).exp();
        }
        base * (y * (delta / u).ln_1p()).exp()
    } else {
        (1.0 + x).powf(y)
    }
}

/// Evaluates the natural log of the binomial PMF $\ln P(X = x)$ for $X \sim \operatorname{Binomial}(n, p)$ with $p + q = 1$.
#[inline]
pub fn log_binomial_pmf(x: f64, n: f64, p: f64, q: f64) -> f64 {
    binomial_pmf(x, n, p, q, true)
}

/// Evaluates the binomial PMF (or log-PMF when `log_p` is `true`) via Loader's (2000, Eq. 6) saddlepoint formula.
#[inline]
pub fn binomial_pmf(x: f64, n: f64, p: f64, q: f64, log_p: bool) -> f64 {
    if x.is_nan() || n.is_nan() || p.is_nan() || q.is_nan() || n < 0.0 || p < 0.0 || q < 0.0 {
        return f64::NAN;
    }
    if x < 0.0 || x > n {
        return prob_zero(log_p);
    }
    if n == 0.0 {
        return prob_one(log_p);
    }
    if p == 0.0 {
        return if x == 0.0 {
            prob_one(log_p)
        } else {
            prob_zero(log_p)
        };
    }
    if q == 0.0 {
        return if x == n {
            prob_one(log_p)
        } else {
            prob_zero(log_p)
        };
    }
    if x == 0.0 {
        return if p <= q {
            if log_p {
                n * (-p).ln_1p()
            } else {
                one_plus_pow(-p, n)
            }
        } else if log_p {
            n * q.ln()
        } else {
            q.powf(n)
        };
    }
    if x == n {
        return if p > q {
            if log_p {
                n * (-q).ln_1p()
            } else {
                one_plus_pow(-q, n)
            }
        } else if log_p {
            n * p.ln()
        } else {
            p.powf(n)
        };
    }
    let n_minus_x = n - x;
    let lc = stirling_error(n)
        - stirling_error(x)
        - stirling_error(n_minus_x)
        - binomial_deviance(x, n * p)
        - binomial_deviance(n_minus_x, n * q);
    let lf = LN_2PI + x.ln() + (-x / n).ln_1p();
    let lp = lc - 0.5 * lf;
    if log_p {
        lp
    } else {
        lp.exp()
    }
}

/// Evaluates the hypergeometric PMF (or log-PMF when `log_p` is `true`) via Loader's (2000, §5) three-binomial composition.
#[inline]
pub fn hypergeom_pmf(x: f64, r: f64, b: f64, n: f64, log_p: bool) -> f64 {
    if x.is_nan() || r.is_nan() || b.is_nan() || n.is_nan() {
        return f64::NAN;
    }
    if r < 0.0 || b < 0.0 || n < 0.0 || !(r + b).is_finite() || !n.is_finite() {
        return f64::NAN;
    }
    let r_int = r.round();
    let b_int = b.round();
    let n_int = n.round();
    if (r - r_int).abs() > 1e-7
        || (b - b_int).abs() > 1e-7
        || (n - n_int).abs() > 1e-7
        || n_int > r_int + b_int
    {
        return f64::NAN;
    }
    if !x.is_finite() {
        return prob_zero(log_p);
    }
    let x_int = x.round();
    if (x - x_int).abs() > 1e-7 {
        return prob_zero(log_p);
    }
    if x_int < 0.0 || x_int > r_int || x_int > n_int || n_int - x_int > b_int {
        return prob_zero(log_p);
    }
    if n_int == 0.0 {
        return prob_one(log_p);
    }
    let total = r_int + b_int;
    let p = n_int / total;
    let q = (total - n_int) / total;
    if log_p {
        let p1 = log_binomial_pmf(x_int, r_int, p, q);
        let p2 = log_binomial_pmf(n_int - x_int, b_int, p, q);
        let p3 = log_binomial_pmf(n_int, total, p, q);
        return p1 + p2 - p3;
    }
    let b1 = binomial_pmf(x_int, r_int, p, q, false);
    let b2 = binomial_pmf(n_int - x_int, b_int, p, q, false);
    let b3 = binomial_pmf(n_int, total, p, q, false);
    if b3 <= 0.0 || b1 == 0.0 || b2 == 0.0 {
        let lp = log_binomial_pmf(x_int, r_int, p, q)
            + log_binomial_pmf(n_int - x_int, b_int, p, q)
            - log_binomial_pmf(n_int, total, p, q);
        return lp.exp().clamp(0.0, 1.0);
    }
    ((b1 * b2) / b3).clamp(0.0, 1.0)
}

/// Evaluates the hypergeometric CDF ($P(X \le x)$ when `lower_tail` is `true`, $P(X > x)$ otherwise)
/// via upper-tail reflection and Abramowitz & Stegun (`§26.3.29`) downward ratio summation.
pub fn hypergeom_cdf(x: f64, nr: f64, nb: f64, n: f64, mut lower_tail: bool) -> f64 {
    if x.is_nan() || nr.is_nan() || nb.is_nan() || n.is_nan() {
        return f64::NAN;
    }
    let mut x = (x + 1e-7).floor();
    let mut nr = nr.round();
    let mut nb = nb.round();
    let n = n.round();
    if nr < 0.0 || nb < 0.0 || !(nr + nb).is_finite() || n < 0.0 || n > nr + nb {
        return f64::NAN;
    }
    if x < 0.0 || x < n - nb {
        return if lower_tail { 0.0 } else { 1.0 };
    }
    if x >= nr || x >= n {
        return if lower_tail { 1.0 } else { 0.0 };
    }
    let reflect = {
        let lhs = x * (nr + nb);
        let rhs = n * nr;
        if lhs.is_finite() && rhs.is_finite() {
            lhs > rhs
        } else {
            x / n > nr / (nr + nb)
        }
    };
    if reflect {
        core::mem::swap(&mut nr, &mut nb);
        x = n - x - 1.0;
        lower_tail = !lower_tail;
        if x < 0.0 || x < n - nb {
            return if lower_tail { 0.0 } else { 1.0 };
        }
        if x >= nr || x >= n {
            return if lower_tail { 1.0 } else { 0.0 };
        }
    }
    let pd = hypergeom_pmf(x, nr, nb, n, false);
    if pd == 0.0 {
        return if lower_tail { 0.0 } else { 1.0 };
    }
    let mut term = 1.0_f64;
    let mut tail_rel = 0.0_f64;
    let mut k = x;
    while k > 0.0 && k - 1.0 < k {
        let numer = k * (nb - n + k);
        if numer <= 0.0 {
            break;
        }
        let denom = (n + 1.0 - k) * (nr + 1.0 - k);
        term *= numer / denom;
        tail_rel += term;
        if term < f64::EPSILON * (1.0 + tail_rel) {
            break;
        }
        k -= 1.0;
    }
    let lower_prob = (pd * (1.0 + tail_rel)).clamp(0.0, 1.0);
    if lower_tail {
        lower_prob
    } else {
        (0.5 - lower_prob + 0.5).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stirling_error_domain_nan_negative_and_infinity() {
        assert!(stirling_error(f64::NAN).is_nan());
        assert!(stirling_error(-1e-15).is_nan());
        assert!(stirling_error(-1.0).is_nan());
        assert!(stirling_error(f64::NEG_INFINITY).is_nan());
        assert_eq!(stirling_error(f64::INFINITY), 0.0);
        assert_eq!(stirling_error(0.0), 0.0);
        assert_eq!(stirling_error(-0.0), 0.0);
    }

    #[test]
    fn stirling_error_half_integers_and_large_n_asymptotic() {
        let half_ln_2pi = 0.5 * LN_2PI;
        let mut ln_fact = 0.0_f64;
        for m in 1..=15 {
            let mf = f64::from(m);
            ln_fact += mf.ln();
            let exact = ln_fact - (mf + 0.5) * mf.ln() + mf - half_ln_2pi;
            let got = stirling_error(mf);
            assert!(
                (got - exact).abs() < 1e-14,
                "stirling_error({mf}): got {got}, expected {exact}"
            );
        }
        // For m = 15..=500, verify the cancellation-free exact step identity
        // delta(m) - delta(m + 1) = (m + 0.5) * ln(1 + 1/m) - 1.
        for m in 15..=500 {
            let mf = f64::from(m);
            let lhs = stirling_error(mf) - stirling_error(mf + 1.0);
            let rhs = (mf + 0.5) * (1.0 / mf).ln_1p() - 1.0;
            assert!(
                (lhs - rhs).abs() < 1e-15,
                "step identity mismatch at m={m}: lhs={lhs}, rhs={rhs}"
            );
        }

        let exact_half =
            (0.5 * std::f64::consts::PI.sqrt()).ln() - 0.5_f64.ln() + 0.5 - half_ln_2pi;
        assert!((stirling_error(0.5) - exact_half).abs() < 1e-15);

        for &big in &[2.0e7 + 1.0, 1.0e9, 1.0e15] {
            assert_eq!(
                stirling_error(big).to_bits(),
                ((1.0 / 12.0) / big).to_bits()
            );
        }
        // At n = 2.0e7, the 8-term Euler-Maclaurin series matches (1/12)/n within <= 1 ULP.
        assert!(
            stirling_error(2.0e7)
                .to_bits()
                .abs_diff(((1.0_f64 / 12.0) / 2.0e7).to_bits())
                <= 1
        );
    }

    #[test]
    fn stirling_error_non_half_integer_shift16_recurrence_and_continuity() {
        for &z in &[
            0.05,
            0.25,
            0.73,
            1.25,
            std::f64::consts::PI,
            7.8,
            12.4,
            14.3,
            14.99,
        ] {
            let lhs = stirling_error(z) - stirling_error(z + 1.0);
            let rhs = (z + 0.5) * (1.0 / z).ln_1p() - 1.0;
            assert!(
                (lhs - rhs).abs() < 1e-13,
                "recurrence mismatch at z={z}: lhs={lhs}, rhs={rhs}"
            );
        }

        for k in 1..=30 {
            let z = f64::from(k) * 0.5;
            let mid = stirling_error(z);
            let minus = stirling_error(z - 1e-9);
            let plus = stirling_error(z + 1e-9);
            assert!(
                (minus - mid).abs() < 1e-8 && (plus - mid).abs() < 1e-8,
                "discontinuity near z={z}: minus={minus}, mid={mid}, plus={plus}"
            );
        }
    }

    #[test]
    fn binomial_deviance_domain_series_cancellation_and_extremes() {
        assert!(binomial_deviance(f64::NAN, 1.0).is_nan());
        assert!(binomial_deviance(1.0, f64::NAN).is_nan());
        assert!(binomial_deviance(f64::INFINITY, 1.0).is_nan());
        assert!(binomial_deviance(1.0, f64::INFINITY).is_nan());
        assert!(binomial_deviance(-0.1, 1.0).is_nan());
        assert!(binomial_deviance(1.0, 0.0).is_nan());
        assert!(binomial_deviance(1.0, -2.0).is_nan());
        assert_eq!(binomial_deviance(0.0, 3.75).to_bits(), 3.75_f64.to_bits());
        assert_eq!(binomial_deviance(4.5, 4.5), 0.0);

        let mu = 10.0;
        let x_hi_in = (11.0 / 9.0) * mu * (1.0 - 1e-12);
        let x_hi_out = (11.0 / 9.0) * mu * (1.0 + 1e-12);
        assert!((binomial_deviance(x_hi_in, mu) - binomial_deviance(x_hi_out, mu)).abs() < 1e-11);
        let x_lo_in = (9.0 / 11.0) * mu * (1.0 + 1e-12);
        let x_lo_out = (9.0 / 11.0) * mu * (1.0 - 1e-12);
        assert!((binomial_deviance(x_lo_in, mu) - binomial_deviance(x_lo_out, mu)).abs() < 1e-11);

        let x_e = std::f64::consts::E * (1.0 + 1e-8);
        let got_e = binomial_deviance(x_e, 1.0);
        let expected_e = x_e * (x_e.ln() - 1.0) + 1.0;
        assert!((got_e - expected_e).abs() < 1e-15);

        let big_in = binomial_deviance(1.5e308, 1.4e308);
        assert!(big_in.is_finite() && big_in > 0.0);
        let big_out = binomial_deviance(1.7e308, 1.0e308);
        assert!(big_out.is_finite() && big_out > 0.0);

        let sub_a = f64::from_bits(4);
        let sub_b = f64::from_bits(5);
        assert!(binomial_deviance(sub_a, sub_b).is_finite());
    }

    #[test]
    fn one_plus_pow_and_binomial_pmf_boundaries_and_identities() {
        assert!(one_plus_pow(f64::NAN, 2.0).is_nan());
        assert!(one_plus_pow(0.1, f64::NAN).is_nan());
        assert_eq!(one_plus_pow(0.25, 0.0), 1.0);
        assert_eq!(one_plus_pow(0.25, 1.0), 1.25);
        assert!((one_plus_pow(0.25, 2.0) - 1.5625).abs() < 1e-15);
        assert!((one_plus_pow(0.25, 3.0) - 1.953125).abs() < 1e-15);
        assert!((one_plus_pow(0.25, 4.0) - 2.44140625).abs() < 1e-15);
        assert_eq!(one_plus_pow(-1.0, 3.0), 0.0);
        let tiny = one_plus_pow(-1e-15, 1e12);
        assert!((tiny - (-1e-3_f64).exp()).abs() < 1e-14);

        assert_eq!(binomial_pmf(0.0, 10.0, 0.0, 1.0, false), 1.0);
        assert_eq!(binomial_pmf(0.0, 10.0, 0.0, 1.0, true), 0.0);
        assert_eq!(binomial_pmf(1.0, 10.0, 0.0, 1.0, false), 0.0);
        assert_eq!(binomial_pmf(1.0, 10.0, 0.0, 1.0, true), f64::NEG_INFINITY);
        assert_eq!(binomial_pmf(10.0, 10.0, 1.0, 0.0, false), 1.0);
        assert_eq!(binomial_pmf(10.0, 10.0, 1.0, 0.0, true), 0.0);
        assert_eq!(binomial_pmf(9.0, 10.0, 1.0, 0.0, false), 0.0);
        assert_eq!(binomial_pmf(9.0, 10.0, 1.0, 0.0, true), f64::NEG_INFINITY);
        assert_eq!(binomial_pmf(-1.0, 10.0, 0.3, 0.7, false), 0.0);
        assert_eq!(binomial_pmf(11.0, 10.0, 0.3, 0.7, true), f64::NEG_INFINITY);
        assert_eq!(binomial_pmf(0.0, 0.0, 0.3, 0.7, false), 1.0);

        for &(n, p) in &[(10, 0.01), (25, 0.25), (50, 0.5), (100, 0.9)] {
            let q = 1.0 - p;
            let nf = f64::from(n);
            let sum: f64 = (0..=n)
                .map(|k| binomial_pmf(f64::from(k), nf, p, q, false))
                .sum();
            assert!(
                (sum - 1.0).abs() < 1e-13,
                "binomial sum={sum} for n={n}, p={p}"
            );
            for k in 0..=n {
                let kf = f64::from(k);
                assert_eq!(
                    log_binomial_pmf(kf, nf, p, q).to_bits(),
                    binomial_pmf(kf, nf, p, q, true).to_bits()
                );
            }
        }
    }

    #[test]
    fn hypergeom_pmf_and_cdf_reference_identities_and_2pow53_termination() {
        assert!(hypergeom_pmf(f64::NAN, 5.0, 5.0, 5.0, false).is_nan());
        assert!(hypergeom_cdf(f64::NAN, 5.0, 5.0, 5.0, true).is_nan());
        assert!(hypergeom_cdf(2.0, -1.0, 5.0, 3.0, true).is_nan());
        assert!(hypergeom_cdf(2.0, 5.0, -1.0, 3.0, true).is_nan());
        assert!(hypergeom_cdf(2.0, 5.0, 5.0, -1.0, true).is_nan());
        assert!(hypergeom_cdf(2.0, 5.0, 5.0, 11.0, true).is_nan());
        assert!(hypergeom_cdf(2.0, f64::INFINITY, 5.0, 3.0, true).is_nan());

        assert_eq!(hypergeom_cdf(f64::INFINITY, 0.0, 0.0, 0.0, true), 1.0);
        assert_eq!(hypergeom_cdf(f64::INFINITY, 0.0, 0.0, 0.0, false), 0.0);
        assert_eq!(hypergeom_cdf(f64::NEG_INFINITY, 5.0, 5.0, 3.0, true), 0.0);
        assert_eq!(hypergeom_cdf(f64::NEG_INFINITY, 5.0, 5.0, 3.0, false), 1.0);

        let exact_counts = [6.0, 60.0, 120.0, 60.0, 6.0];
        let mut running = 0.0;
        for (x, &num) in exact_counts.iter().enumerate() {
            let xf = x as f64;
            let expected_pmf = num / 252.0;
            running += expected_pmf;
            let got_pmf = hypergeom_pmf(xf, 4.0, 6.0, 5.0, false);
            assert!((got_pmf - expected_pmf).abs() < 1e-14);
            let sym_pmf = hypergeom_pmf(5.0 - xf, 6.0, 4.0, 5.0, false);
            assert!((got_pmf - sym_pmf).abs() < 1e-15);

            let cdf_lo = hypergeom_cdf(xf, 4.0, 6.0, 5.0, true);
            let cdf_hi = hypergeom_cdf(xf, 4.0, 6.0, 5.0, false);
            assert!((cdf_lo - running).abs() < 1e-14);
            assert!((cdf_lo + cdf_hi - 1.0).abs() < 1e-14);
        }

        let t0 = std::time::Instant::now();
        let two53 = 9_007_199_254_740_992.0_f64;
        assert_eq!(hypergeom_cdf(two53, 100.0, 100.0, 50.0, true), 1.0);
        assert_eq!(hypergeom_cdf(two53, 100.0, 100.0, 50.0, false), 0.0);
        let c1 = hypergeom_cdf(
            two53 + 4.0,
            2.0_f64.powi(55),
            2.0_f64.powi(55),
            2.0_f64.powi(55),
            true,
        );
        let c2 = hypergeom_cdf(
            2.0_f64.powi(54),
            2.0_f64.powi(60),
            2.0_f64.powi(60),
            2.0_f64.powi(55),
            true,
        );
        assert!(c1.is_finite() && (0.0..=1.0).contains(&c1));
        assert!(c2.is_finite() && (0.0..=1.0).contains(&c2));
        assert!(
            t0.elapsed() < std::time::Duration::from_millis(200),
            "x >= 2^53 termination took too long: {:?}",
            t0.elapsed()
        );
    }
}
