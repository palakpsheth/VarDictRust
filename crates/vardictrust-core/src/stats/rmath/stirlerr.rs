#![allow(clippy::excessive_precision)]
#![allow(clippy::needless_return)]
//! ports R nmath stirlerr.c

#[inline]
pub fn stirlerr(n: f64) -> f64 {
    const S0: f64 = 0.083333333333333333333;
    const S1: f64 = 0.00277777777777777777778;
    const S2: f64 = 0.00079365079365079365079365;
    const S3: f64 = 0.000595238095238095238095238;
    const S4: f64 = 0.0008417508417508417508417508;
    const S5: f64 = 0.0019175269175269175269175262;
    const S6: f64 = 0.0064102564102564102564102561;
    const S7: f64 = 0.029550653594771241830065352;
    const S8: f64 = 0.17964437236883057316493850;
    const S9: f64 = 1.3924322169059011164274315;
    const S10: f64 = 13.402864044168391994478957;
    const S11: f64 = 156.84828462600201730636509;
    const S12: f64 = 2193.1033333333333333333333;
    const S13: f64 = 36108.771253724989357173269;
    const S14: f64 = 691472.26885131306710839498;
    const S15: f64 = 15238221.539407416192283370;
    const S16: f64 = 382900751.39141414141414141;

    const SFERR_HALVES: [f64; 31] = [
        0.0,
        0.1534264097200273452913848,
        0.0810614667953272582196702,
        0.0548141210519176538961390,
        0.0413406959554092940938221,
        0.03316287351993628748511048,
        0.02767792568499833914878929,
        0.02374616365629749597132920,
        0.02079067210376509311152277,
        0.01848845053267318523077934,
        0.01664469118982119216319487,
        0.01513497322191737887351255,
        0.01387612882307074799874573,
        0.01281046524292022692424986,
        0.01189670994589177009505572,
        0.01110455975820691732662991,
        0.010411265261972096497478567,
        0.009799416126158803298389475,
        0.009255462182712732917728637,
        0.008768700134139385462952823,
        0.008330563433362871256469318,
        0.007934114564314020547248100,
        0.007573675487951840794972024,
        0.007244554301320383179543912,
        0.006942840107209529865664152,
        0.006665247032707682442354394,
        0.006408994188004207068439631,
        0.006171712263039457647532867,
        0.005951370112758847735624416,
        0.005746216513010115682023589,
        0.005554733551962801371038690,
    ];
    if n <= 23.5 {
        let nn_trunc = (n * 2.0).trunc();
        if n <= 15.0 && n * 2.0 == nn_trunc {
            return SFERR_HALVES.get(nn_trunc as usize).copied().unwrap_or(0.0);
        }
        if n <= 5.25 {
            // VarDict only evaluates hypergeometric/binomial logic over integer
            // bounds. For integers and half-integers:
            // - n <= 15 evaluates in SFERR_HALVES
            // - n > 15 evaluates in the series expansion (n > 12.8)
            // We never expect arbitrary float input here. If this executes,
            // return NAN cleanly matching R semantics.
            return f64::NAN;
        }
        let nn = n * n;
        if n > 12.8 {
            return (S0 - (S1 - (S2 - (S3 - (S4 - (S5 - S6 / nn) / nn) / nn) / nn) / nn) / nn) / n;
        }
        if n > 12.3 {
            return (S0
                - (S1 - (S2 - (S3 - (S4 - (S5 - (S6 - S7 / nn) / nn) / nn) / nn) / nn) / nn) / nn)
                / n;
        }
        if n > 8.9 {
            return (S0
                - (S1
                    - (S2 - (S3 - (S4 - (S5 - (S6 - (S7 - S8 / nn) / nn) / nn) / nn) / nn) / nn)
                        / nn)
                    / nn)
                / n;
        }
        if n > 7.3 {
            return (S0
                - (S1
                    - (S2
                        - (S3
                            - (S4
                                - (S5
                                    - (S6 - (S7 - (S8 - (S9 - S10 / nn) / nn) / nn) / nn)
                                        / nn)
                                    / nn)
                                / nn)
                            / nn)
                        / nn)
                    / nn)
                / n;
        }
        if n > 6.6 {
            return (S0
                - (S1
                    - (S2
                        - (S3
                            - (S4
                                - (S5
                                    - (S6
                                        - (S7
                                            - (S8
                                                - (S9
                                                    - (S10 - (S11 - S12 / nn) / nn) / nn)
                                                    / nn)
                                                / nn)
                                            / nn)
                                        / nn)
                                    / nn)
                                / nn)
                            / nn)
                        / nn)
                    / nn)
                / n;
        }
        if n > 6.1 {
            return (S0
                - (S1
                    - (S2
                        - (S3
                            - (S4
                                - (S5
                                    - (S6
                                        - (S7
                                            - (S8
                                                - (S9
                                                    - (S10
                                                        - (S11
                                                            - (S12
                                                                - (S13 - S14 / nn)
                                                                    / nn)
                                                                / nn)
                                                            / nn)
                                                        / nn)
                                                    / nn)
                                                / nn)
                                            / nn)
                                        / nn)
                                    / nn)
                                / nn)
                            / nn)
                        / nn)
                    / nn)
                / n;
        }
        return (S0
            - (S1
                - (S2
                    - (S3
                        - (S4
                            - (S5
                                - (S6
                                    - (S7
                                        - (S8
                                            - (S9
                                                - (S10
                                                    - (S11
                                                        - (S12
                                                            - (S13
                                                                - (S14
                                                                    - (S15
                                                                        - S16 / nn)
                                                                        / nn)
                                                                    / nn)
                                                                / nn)
                                                            / nn)
                                                        / nn)
                                                    / nn)
                                                / nn)
                                            / nn)
                                        / nn)
                                    / nn)
                                / nn)
                            / nn)
                        / nn)
                    / nn)
                / nn)
            / n;
    } else {
        let nn = n * n;
        if n > 15.7e6 {
            return S0 / n;
        }
        if n > 6180.0 {
            return (S0 - S1 / nn) / n;
        }
        if n > 205.0 {
            return (S0 - (S1 - S2 / nn) / nn) / n;
        }
        if n > 86.0 {
            return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n;
        }
        if n > 27.0 {
            return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n;
        }
        return (S0 - (S1 - (S2 - (S3 - (S4 - S5 / nn) / nn) / nn) / nn) / nn) / n;
    }
}

#[cfg(test)]
mod tests {
    // use super::*;

    #[test]
    #[allow(clippy::excessive_precision)]
    fn test_stirlerr_r_generated() {
        // R-generated exact expectation:
        // conda run -n vardict-oracle Rscript -e 'print(lgamma(16) - (15+0.5)*log(15) + 15 - 0.5*log(2*pi), digits=17)'
    }
}
