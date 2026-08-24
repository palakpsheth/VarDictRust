//! Brent root finder algorithm.
//!
//! Ports `UnirootZeroIn.java`.

/// Finds a root of the function `f` in the interval `[ax, bx]` using Brent's method.
///
/// Ports `UnirootZeroIn.zeroinC`.
pub fn zeroin_c<F>(ax: f64, bx: f64, f: F, tol: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    let epsilon = f64::EPSILON;
    let mut a = ax;
    let mut b = bx;
    let mut fa = f(a);
    let mut fb = f(b);
    let mut c = a;
    let mut fc = fa;

    loop {
        let prev_step = b - a;

        // Swap data for b to be the best approximation
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol_act = 2.0 * epsilon * b.abs() + tol / 2.0;
        let mut new_step = (c - b) / 2.0;

        // Acceptable approx. is found
        if new_step.abs() <= tol_act || fb == 0.0 {
            return b;
        }

        // Decide if the interpolation can be tried
        if prev_step.abs() >= tol_act && fa.abs() > fb.abs() {
            let mut p;
            let mut q;
            let cb = c - b;

            if (a - c).abs() < f64::MIN_POSITIVE {
                // Linear interpolation
                let t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                // Quadratic inverse interpolation
                q = fa / fc;
                let t1 = fb / fc;
                let t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }

            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }

            if p < (0.75 * cb * q - (tol_act * q).abs() / 2.0) && p < (prev_step * q / 2.0).abs() {
                new_step = p / q;
            }
        }

        // Adjust the step to be not less than tolerance
        if new_step.abs() < tol_act {
            if new_step > 0.0 {
                new_step = tol_act;
            } else {
                new_step = -tol_act;
            }
        }

        a = b;
        fa = fb;
        b += new_step;
        fb = f(b);

        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zeroin_simple() {
        let f = |x: f64| x * x - 4.0;
        let root = zeroin_c(0.0, 3.0, f, 1e-15);
        assert!((root - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_zeroin_trig() {
        let f = |x: f64| x.sin();
        let root = zeroin_c(3.0, 4.0, f, 1e-15);
        assert!((root - std::f64::consts::PI).abs() < 1e-12);
    }
}
