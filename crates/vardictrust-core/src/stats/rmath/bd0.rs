//! ports R nmath bd0.c
#[inline]
pub fn bd0(x: f64, np: f64) -> f64 {
    if !x.is_finite() || !np.is_finite() || np == 0.0 {
        return f64::NAN;
    }
    if (x - np).abs() < 0.1 * (x + np) {
        let d = x - np;
        let mut v = d / (x + np);
        if d != 0.0 && v == 0.0 {
            let x_ = x * 0.25;
            let n_ = np * 0.25;
            v = (x_ - n_) / (x_ + n_);
        }
        let mut s = (d * 0.5) * v;
        if (s * 2.0).abs() < f64::MIN_POSITIVE {
            return s * 2.0;
        }
        let mut ej = x * v;
        v *= v;
        for j in 1..1000 {
            ej *= v;
            let s_ = s;
            s += ej / ((j << 1) + 1) as f64;
            if s == s_ {
                return s * 2.0;
            }
        }
        return s * 2.0;
    }
    let lg_x_n = if (x / np).is_finite() {
        (x / np).ln()
    } else {
        x.ln() - np.ln()
    };
    if x > np {
        x * (lg_x_n - 1.0) + np
    } else {
        x * lg_x_n + np - x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bd0_branches() {
        // Asymptotic branch
        // conda run -n vardict-oracle Rscript -e 'print(5 * log(5/2) + 2 - 5, digits=17)'
        // 5 * log(2.5) - 3 = 1.5814536593707753
        assert!((bd0(5.0, 2.0) - 1.5814536593707753).abs() < 1e-15);

        // Series branch: |x-np| < 0.1*(x+np)
        // x = 1.05, np = 1.0
        assert!(bd0(1.05, 1.0) > 0.0);
    }
}
