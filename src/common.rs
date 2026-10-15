use std::f64::consts::PI;

use scirs2::{
    integrate::quad::{QuadOptions, quad},
    special::elliptic_ext::{ellipe_agm, ellipk_agm},
};

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

// let (a_star, b_star, delta_star, e, kappa) = hertz_coff(f_rho)
#[allow(non_snake_case)]
pub fn hertz_coefficients(f_rho: f64) -> anyhow::Result<(f64, f64, f64, f64, f64)> {
    let lo = 1e-5;
    let hi = 1.0 - 1e-5;

    let helper = |e: f64| -> f64 {
        let m = e.powf(2.0);
        let K = ellipk_agm(m).unwrap();
        let E = ellipe_agm(m).unwrap();
        1.0 - 2.0 * (1.0 - m) / m * (K - E) / E - f_rho
    };
    let mut convergency = roots::SimpleConvergency {
        eps: 1e-15,
        max_iter: 100,
    };
    let e = roots::find_root_brent(lo, hi, helper, &mut convergency)?;
    let m = e * e;
    let K = ellipk_agm(m)?;
    let E = ellipe_agm(m)?;
    let a_star = (2.0 * E / (PI * (1.0 - m))).powf(1.0 / 3.0);
    let b_star = (2.0 * (1.0 - m).sqrt() * E / PI).powf(1.0 / 3.0);
    let delta_star = 2.0 * K / PI * (PI * (1.0 - m) / 2.0 / E).powf(1.0 / 3.0);

    Ok((a_star, b_star, delta_star, e, a_star / b_star))
}

#[cfg(test)]
mod tests {
    use super::*;
    macro_rules! assert_approx_eq {
        ($a:expr, $b:expr) => {
            assert_approx_eq!($a, $b, 1e-3);
        };
        ($a:expr, $b:expr, $eps:expr) => {{
            let a = $a;
            let b = $b;
            let eps:f64 = $eps;
            let diff = (a - b).abs();
            let threshold = eps.max(eps * a.abs().max(b.abs()));
            assert!(
                diff < threshold,
                "assert_approx_eq failed:\n  a = {a:.15}\n  b = {b:.15}\n  |a-b| = {diff:.15}\n  threshold = {threshold:.15}"
            );
        }};
    }
    #[test]
    fn test_jra_0_5() -> anyhow::Result<()> {
        // roller
        let (jr, ja) = jra(0.5, 1.1)?;
        assert_approx_eq!(0.2453, jr);
        assert_approx_eq!(0.3090, ja);
        // ball
        let (jr, ja) = jra(0.5, 1.5)?;
        assert_approx_eq!(0.2288, jr);
        assert_approx_eq!(0.2782, ja);
        Ok(())
    }

    #[test]
    fn test_hertz() -> anyhow::Result<()> {
        let (a_star, b_star, delta_star, e, _) = hertz_coefficients(0.5)?;
        assert_approx_eq!(a_star, 1.4858);
        assert_approx_eq!(b_star, 0.7169);
        assert_approx_eq!(delta_star, 0.9376);
        assert_approx_eq!(e, 0.8755);
        Ok(())
    }

    #[test]
    #[ignore]
    fn print_hertz_coefficients() -> anyhow::Result<()> {
        let mut f_rho = 0.05;
        let step = 0.05;
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            "f_rho", "a_star", "b_star", "delta_star", "e", "kappa"
        );
        while f_rho < 1.0 {
            let (a_star, b_star, delta_star, e, kappa) = hertz_coefficients(f_rho)?;
            println!(
                "{:.4}\t{:.4}\t{:.4}\t{:.4}\t\t{:.4}\t{:.4}",
                f_rho, a_star, b_star, delta_star, e, kappa
            );
            f_rho += step;
        }
        Ok(())
    }
}
