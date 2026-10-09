#![allow(clippy::excessive_precision)]
//! ports R nmath dbinom.c

#[inline]
pub fn dbinom_raw(x: f64, n: f64, p: f64, q: f64, give_log: bool) -> f64 {
    if p == 0.0 {
        return if x == 0.0 {
            if give_log {
                0.0
            } else {
                1.0
            }
        } else {
            if give_log {
                f64::NEG_INFINITY
            } else {
                0.0
            }
        };
    }
    if q == 0.0 {
        return if x == n {
            if give_log {
                0.0
            } else {
                1.0
            }
        } else {
            if give_log {
                f64::NEG_INFINITY
            } else {
                0.0
            }
        };
    }
    if x == 0.0 {
        if n == 0.0 {
            return if give_log { 0.0 } else { 1.0 };
        }
        if p > q {
            return if give_log { n * q.ln() } else { q.powf(n) };
        } else {
            return if give_log {
                n * (-p).ln_1p()
            } else {
                (1.0 - p).powf(n)
            };
        }
    }
    if x == n {
        if p > q {
            return if give_log {
                n * (-q).ln_1p()
            } else {
                (1.0 - q).powf(n)
            };
        } else {
            return if give_log { n * p.ln() } else { p.powf(n) };
        }
    }
    if x < 0.0 || x > n {
        return if give_log { f64::NEG_INFINITY } else { 0.0 };
    }
    let lc = super::stirlerr::stirlerr(n)
        - super::stirlerr::stirlerr(x)
        - super::stirlerr::stirlerr(n - x)
        - super::bd0::bd0(x, n * p)
        - super::bd0::bd0(n - x, n * q);

    let lf = 1.8378770664093454835606594728112352797 + x.ln() + (-x / n).ln_1p();
    let res = lc - 0.5 * lf;
    if give_log {
        res
    } else {
        res.exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_dbinom_raw_r_generated() {
        // R-generated exact expectation:
        // conda run -n vardict-oracle Rscript -e 'print(dbinom(5, 10, 0.5, log=FALSE), digits=17)'
        // R version 4.6.0
        assert!((dbinom_raw(5.0, 10.0, 0.5, 0.5, false) - 0.24609375000000003).abs() < 1e-15);
        assert!((dbinom_raw(5.0, 10.0, 0.5, 0.5, true) - (-1.4020427180880297)).abs() < 1e-15);
    }
}
