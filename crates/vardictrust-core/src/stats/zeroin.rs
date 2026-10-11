//! Brent root finder algorithm.
//!
//! Ports `UnirootZeroIn.java` and R's `zeroin.c`.

/// R's R_zeroin2
pub fn zeroin_c<F>(
    ax: f64,
    bx: f64,
    mut fa: f64,
    mut fb: f64,
    mut f: F,
    tol_ptr: &mut f64,
    maxit_ptr: &mut i32,
) -> f64
where
    F: FnMut(f64) -> f64,
{
    let mut a = ax;
    let mut b = bx;
    let mut c = a;
    let mut fc = fa;
    let mut maxit = *maxit_ptr + 1;
    let tol = *tol_ptr;

    if fa == 0.0 {
        *tol_ptr = 0.0;
        *maxit_ptr = 0;
        return a;
    }
    if fb == 0.0 {
        *tol_ptr = 0.0;
        *maxit_ptr = 0;
        return b;
    }

    while maxit > 0 {
        maxit -= 1;
        let prev_step = b - a;
        let mut p;
        let mut q;
        let mut new_step;

        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol_act = 2.0 * f64::EPSILON * b.abs() + tol / 2.0;
        new_step = (c - b) / 2.0;

        if new_step.abs() <= tol_act || fb == 0.0 {
            *maxit_ptr -= maxit;
            *tol_ptr = (c - b).abs();
            return b;
        }

        if prev_step.abs() >= tol_act && fa.abs() > fb.abs() {
            let cb = c - b;
            if a == c {
                let t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
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

    *tol_ptr = (c - b).abs();
    *maxit_ptr = -1;
    b
}

pub fn uniroot<F>(lower: f64, upper: f64, mut f: F, tol: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    let f_lower = f(lower);
    let f_upper = f(upper);
    let mut tol_val = tol;
    let mut maxit = 1000;
    zeroin_c(lower, upper, f_lower, f_upper, f, &mut tol_val, &mut maxit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zeroin_simple() {
        let f = |x: f64| x * x - 4.0;
        let root = uniroot(0.0, 3.0, f, 1e-15);
        assert!((root - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_zeroin_trig() {
        let f = |x: f64| x.sin();
        let root = uniroot(3.0, 4.0, f, 1e-15);
        assert!((root - std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_zeroin_chr16() {
        let fisher = crate::stats::FisherExact::new(5, 21, 1, 12);
        let eval = |t: f64| fisher.mnhyper(1.0 / t) - 5.0;
        let eps = f64::EPSILON;
        let tol = eps.powf(0.25);
        let root = uniroot(eps, 1.0, eval, tol);
        assert!((root - 0.35818903155718163).abs() < 1e-12);
    }
}
