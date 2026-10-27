pub mod bd0;
pub mod dbinom_raw;
pub mod dhyper;
pub mod phyper;
pub mod stirlerr;

#[cfg(test)]
mod tests {
    use super::bd0::bd0;
    use super::dbinom_raw::{dbinom_raw, pow1p};
    use super::dhyper::dhyper;
    use super::phyper::phyper;
    use super::stirlerr::stirlerr;

    #[test]
    fn rmath_boundary_cases() {
        // bd0 non-finite / zero np / exact equality / series / asymptotic branches
        assert!(bd0(f64::NAN, 1.0).is_nan());
        assert!(bd0(1.0, f64::INFINITY).is_nan());
        assert!(bd0(1.0, 0.0).is_nan());
        assert_eq!(bd0(5.0, 5.0), 0.0);
        assert!(bd0(10.0, 2.0) > 0.0);
        assert!(bd0(2.0, 10.0) > 0.0);

        // stirlerr half-integer table, non-half-integer <= 5.25 NaN, and polynomial tiers
        assert_eq!(stirlerr(0.0), 0.0);
        assert!(stirlerr(0.5) > 0.15);
        assert!(stirlerr(3.1).is_nan());
        for &n in &[
            5.5, 6.3, 7.0, 8.0, 10.0, 12.5, 16.0, 30.0, 50.0, 100.0, 500.0, 10_000.0, 2e7,
        ] {
            let s = stirlerr(n);
            assert!(s.is_finite() && s > 0.0, "stirlerr({n}) = {s}");
        }

        // pow1p boundary branches
        assert_eq!(pow1p(0.0, f64::NAN), 1.0);
        assert!(pow1p(0.5, f64::NAN).is_nan());
        for y in 0..=4 {
            let expected = (1.25f64).powi(y);
            assert!((pow1p(0.25, y as f64) - expected).abs() < 1e-15);
        }
        assert!((pow1p(1e-12, 2.5) - (1.0 + 2.5e-12)).abs() < 1e-18);
        assert!((pow1p(0.75, 2.5) - 1.75f64.powf(2.5)).abs() < 1e-15);

        // dbinom_raw boundary branches (p=0, q=0, x=0, x=n, out-of-range)
        assert_eq!(dbinom_raw(0.0, 5.0, 0.0, 1.0, false), 1.0);
        assert_eq!(dbinom_raw(0.0, 5.0, 0.0, 1.0, true), 0.0);
        assert_eq!(dbinom_raw(1.0, 5.0, 0.0, 1.0, false), 0.0);
        assert_eq!(dbinom_raw(1.0, 5.0, 0.0, 1.0, true), f64::NEG_INFINITY);
        assert_eq!(dbinom_raw(5.0, 5.0, 1.0, 0.0, false), 1.0);
        assert_eq!(dbinom_raw(4.0, 5.0, 1.0, 0.0, false), 0.0);
        assert_eq!(dbinom_raw(0.0, 0.0, 0.5, 0.5, false), 1.0);
        assert_eq!(dbinom_raw(-1.0, 5.0, 0.5, 0.5, false), 0.0);
        assert_eq!(dbinom_raw(6.0, 5.0, 0.5, 0.5, true), f64::NEG_INFINITY);

        // dhyper boundary branches
        assert_eq!(dhyper(-1.0, 5.0, 5.0, 5.0, false), 0.0);
        assert_eq!(dhyper(6.0, 5.0, 5.0, 5.0, true), f64::NEG_INFINITY);
        assert_eq!(dhyper(0.0, 5.0, 5.0, 0.0, false), 1.0);
        assert_eq!(dhyper(0.0, 5.0, 5.0, 0.0, true), 0.0);

        // phyper NaN, invalid domain, tail-swap, and bounds
        assert!(phyper(f64::NAN, 5.0, 5.0, 5.0, true).is_nan());
        assert!(phyper(1.0, -1.0, 5.0, 5.0, true).is_nan());
        assert!(phyper(1.0, 5.0, 5.0, 12.0, true).is_nan());
        assert_eq!(phyper(-1.0, 5.0, 5.0, 5.0, true), 0.0);
        assert_eq!(phyper(-1.0, 5.0, 5.0, 5.0, false), 1.0);
        assert_eq!(phyper(5.0, 5.0, 5.0, 5.0, true), 1.0);
        assert_eq!(phyper(5.0, 5.0, 5.0, 5.0, false), 0.0);
        // Complement identity across tail swap
        let lower = phyper(4.0, 5.0, 5.0, 5.0, true);
        let upper = phyper(4.0, 5.0, 5.0, 5.0, false);
        assert!((lower + upper - 1.0).abs() < 1e-14);
    }
}
