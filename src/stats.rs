//! Chi-squared survival function via the regularised upper incomplete gamma
//! function. P(chi²_df > x) = Q(df/2, x/2).

#[must_use]
pub fn chi2_sf(x: f64, df: u32) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 1.0;
    }
    if x.is_infinite() {
        return 0.0;
    }
    gamma_q(f64::from(df) / 2.0, x / 2.0)
}

fn gamma_q(a: f64, z: f64) -> f64 {
    if z < a + 1.0 {
        1.0 - gamma_p_series(a, z)
    } else {
        gamma_q_cf(a, z)
    }
}

fn gamma_p_series(a: f64, z: f64) -> f64 {
    if z <= 0.0 {
        return 0.0;
    }
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..200 {
        ap += 1.0;
        del *= z / ap;
        sum += del;
        if del.abs() < sum.abs() * 3e-15 {
            break;
        }
    }
    sum * (-z + a * z.ln() - ln_gamma(a)).exp()
}

fn gamma_q_cf(a: f64, z: f64) -> f64 {
    let fpmin = f64::MIN_POSITIVE / 3e-15;
    let mut b = z + 1.0 - a;
    let mut c = 1.0 / fpmin;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1u32..=200 {
        let an = -(f64::from(i) * (f64::from(i) - a));
        b += 2.0;
        d = an * d + b;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = b + an / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 3e-15 {
            break;
        }
    }
    h * (-z + a * z.ln() - ln_gamma(a)).exp()
}

/// Lanczos approximation (g=7, n=9).
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_403,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_5,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * x).sin().ln()
            - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + G + 0.5;
    let ser: f64 = C[0]
        + C[1..]
            .iter()
            .enumerate()
            .fold(0.0, |acc, (i, &c)| acc + c / (x + i as f64 + 1.0));
    0.5 * (2.0 * std::f64::consts::PI).ln() + ser.ln() + (x + 0.5) * t.ln() - t
}
