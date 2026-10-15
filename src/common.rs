use std::f64::consts::PI;

use anyhow::Ok;
use scirs2::integrate::quad::{QuadOptions, quad};

// 载荷分布积分计算
pub fn jra(epsilon: f64, n: f64) -> anyhow::Result<(f64, f64)> {
    let f_a = |psi: f64| -> f64 {
        let term = 1.0 - 1.0 / (2.0 * epsilon) * (1.0 - psi.cos());
        if term <= 0.0 {
            0.0
        } else {
            term.powf(n) / (2.0 * PI)
        }
    };
    let f_r = |psi: f64| f_a(psi) * psi.cos();
    // 积分选项
    let opts = QuadOptions {
        use_simpson: true,
        ..Default::default()
    };
    let ja = quad(f_a, -PI, PI, Some(opts.clone()))?;
    let jr = quad(f_r, -PI, PI, Some(opts))?;
    Ok((jr.value, ja.value))
}


// 使用和matlab相同算法，scirs2中近似算法误差较大
fn ellipke(m: f64) -> (f64, f64) {
    let tol = 1e-15;
    let max_iter = 100;
    let mut a = 1.0_f64;
    let mut b = (1.0 - m).sqrt();
    let mut c = m.sqrt();
    let mut sum = 0.5 * c * c; // n = 0 项，权重 2^(-1)

    let mut n: i32 = 0;
    while c.abs() > tol && n < max_iter {
        let a_next = 0.5 * (a + b);
        let b_next = (a * b).sqrt();
        let c_next = 0.5 * (a - b);

        a = a_next;
        b = b_next;
        c = c_next;
        n += 1;
        sum += 2f64.powi(n - 1) * c * c;
    }

    let k = PI / (2.0 * a);
    let e = k * (1.0 - sum);
    (k, e)
}

// let (a_star, b_star, delta_star, e, kappa) = hertz_coff(f_rho)
#[allow(non_snake_case)]
pub fn hertz_coff(f_rho: f64) -> anyhow::Result<(f64, f64, f64, f64, f64)> {
    let lo = 1e-5;
    let hi = 1.0 - 1e-5;

    let helper = |e: f64| -> f64 {
        let m = e.powf(2.0);
        let (K, E) = ellipke(m);
        1.0 - 2.0 * (1.0 - m) / m * (K - E) / E - f_rho
    };
    let mut convergency = roots::SimpleConvergency {
        eps: 1e-15,
        max_iter: 100,
    };
    let e = roots::find_root_brent(lo, hi, helper, &mut convergency)?;
    let m = e * e;
    let (K, E) = ellipke(m);

    let a_star = (2.0 * E / (PI * (1.0 - m))).powf(1.0 / 3.0);
    let b_star = (2.0 * (1.0 - m).sqrt() * E / PI).powf(1.0 / 3.0);
    let delta_star = 2.0 * K / PI * (PI * (1.0 - m) / 2.0 / E).powf(1.0 / 3.0);

    Ok((a_star, b_star, delta_star, e, a_star / b_star))
}


#[cfg(test)]
fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

#[cfg(test)]
mod tests {
    use anyhow::Ok;

    use super::*;

    #[test]
    fn test_jra_0_5() -> anyhow::Result<()> {
        let eps = 1e-3;
        // roller
        let (jr, ja) = jra(0.5, 1.1)?;
        assert!(approx_eq(0.2453, jr, eps));
        assert!(approx_eq(0.3090, ja, eps));
        // ball
        let (jr, ja) = jra(0.5, 1.5)?;
        assert!(approx_eq(0.2288, jr, eps));
        assert!(approx_eq(0.2782, ja, eps));
        Ok(())
    }

    #[test]
    fn test_hertz() -> anyhow::Result<()> {
        let eps = 1e-3;
        let (a_star, b_star, delta_star,  e, _ )= hertz_coff(0.5)?;
        assert!(approx_eq(a_star, 1.4858, eps));
        assert!(approx_eq(b_star, 0.7169, eps));
        assert!(approx_eq(delta_star, 0.9376, eps));
        assert!(approx_eq(e, 0.8755, eps));
        Ok(())
    }

    #[test]
    fn test_ellipke() {
        let eps = 1e-3;
        let (k, e) = ellipke(0.9);
        assert!(approx_eq(k, 2.5781, eps));
        assert!(approx_eq(e, 1.1048, eps))
    }


    #[test]
    #[ignore]
    fn print_hertz_coefficients()-> anyhow::Result<()>{
        let mut f_rho = 0.05;
        let step = 0.05;
            println!("{}\t{}\t{}\t{}\t{}\t{}", "f_rho", "a_star", "b_star", "delta_star", "e", "kappa");
        while f_rho < 1.0 {
            let (a_star , b_star, delta_star, e,  kappa) = hertz_coff(f_rho)?;
            println!("{:.4}\t{:.4}\t{:.4}\t{:.4}\t\t{:.4}\t{:.4}", f_rho,a_star, b_star, delta_star, e, kappa);
            f_rho += step;
        }
        Ok(())
    }
}
