//! ports R nmath phyper.c

use super::dhyper::dhyper;

#[inline]
fn pdhyper(mut x: f64, nr: f64, nb: f64, n: f64) -> f64 {
    let mut sum = 0.0;
    let mut term = 1.0;

    while x > 0.0 && term >= f64::EPSILON * sum {
        term *= (x * (nb - n + x)) / ((n + 1.0 - x) * (nr + 1.0 - x));
        sum += term;
        x -= 1.0;
    }

    1.0 + sum
}

/// Distribution function of the hypergeometric distribution.
///
/// Ports `Rf_phyper` from R's `nmath/phyper.c`.
#[inline]
pub fn phyper(mut x: f64, mut nr: f64, mut nb: f64, mut n: f64, mut lower_tail: bool) -> f64 {
    if x.is_nan() || nr.is_nan() || nb.is_nan() || n.is_nan() {
        return x + nr + nb + n;
    }

    x = (x + 1e-7).floor();
    nr = nr.round();
    nb = nb.round();
    n = n.round();

    if nr < 0.0 || nb < 0.0 || !(nr + nb).is_finite() || n < 0.0 || n > nr + nb {
        return f64::NAN;
    }

    if x * (nr + nb) > n * nr {
        /* Swap tails. */
        std::mem::swap(&mut nb, &mut nr);
        x = n - x - 1.0;
        lower_tail = !lower_tail;
    }

    if x < 0.0 || x < n - nb {
        return if lower_tail { 0.0 } else { 1.0 };
    }
    if x >= nr || x >= n {
        return if lower_tail { 1.0 } else { 0.0 };
    }

    let d = dhyper(x, nr, nb, n, false);
    if d == 0.0 {
        return if lower_tail { 0.0 } else { 1.0 };
    }
    let pd = pdhyper(x, nr, nb, n);
    let val = d * pd;

    if lower_tail {
        val
    } else {
        1.0 - val
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_phyper_r_generated() {
        let v = phyper(86.0, 87.0, 105.0, 191.0, true);
        assert!((v - 0.45312500000000189).abs() < 1e-14);
        assert!((phyper(85.0, 87.0, 105.0, 191.0, false) - 1.0).abs() < 1e-15);
        let v2 = phyper(0.0, 42.0, 87.0, 2.0, true);
        assert!((v2 - 0.45312500000000011).abs() < 1e-15);
    }
}
