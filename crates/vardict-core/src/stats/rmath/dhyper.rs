//! ports R nmath dhyper.c
#[inline]
pub fn dhyper(x: f64, r: f64, b: f64, n: f64, give_log: bool) -> f64 {
    if x < 0.0 {
        return if give_log { f64::NEG_INFINITY } else { 0.0 };
    }
    if n < x || r < x || n - x > b {
        return if give_log { f64::NEG_INFINITY } else { 0.0 };
    }
    if n == 0.0 {
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
    let p = n / (r + b);
    let q = (r + b - n) / (r + b);
    let p1 = super::dbinom_raw::dbinom_raw(x, r, p, q, give_log);
    let p2 = super::dbinom_raw::dbinom_raw(n - x, b, p, q, give_log);
    let p3 = super::dbinom_raw::dbinom_raw(n, r + b, p, q, give_log);
    if give_log {
        p1 + p2 - p3
    } else {
        p1 * p2 / p3
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_dhyper_r_generated() {
        // R-generated exact expectation:
        // conda run -n vardict-oracle Rscript -e 'print(dhyper(5, 10, 10, 10, log=FALSE), digits=17)'
        // R version 4.6.0
        assert!((dhyper(5.0, 10.0, 10.0, 10.0, false) - 0.34371820130334069).abs() < 1e-15);
        assert!((dhyper(5.0, 10.0, 10.0, 10.0, true) - (-1.0679331395796077)).abs() < 1e-15);
    }
}
