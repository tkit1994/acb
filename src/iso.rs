use anyhow::ensure;
use scirs2::special::elliptic_ext::ellipe_agm;

use crate::{
    common::hertz_coefficients,
    core::{calc_sum_rho_f_rho_i, calc_sum_rho_f_rho_o},
};

// ISO 281 径向额定中载荷系数计算
pub fn fcr(a1: f64, lambda: f64, fi: f64, fe: f64, gamma: f64) -> anyhow::Result<f64> {
    // 1. 使用 anyhow 的 ensure! 宏进行严格的物理边界检查
    ensure!(
        (2.0 * fi - 1.0).abs() >= f64::EPSILON,
        "计算错误: 分母 (2 * fi - 1) 接近零，导致除零异常 (fi = {})",
        fi
    );
    ensure!(
        fe.abs() >= f64::EPSILON,
        "计算错误: 外圈曲率系数 fe 不能为零"
    );
    ensure!(
        (0.0..1.0).contains(&gamma),
        "计算错误: 几何比 gamma 超出有效范围 [0, 1) (实际值: {})",
        gamma
    );

    // 2. 基础常数与线性组合项
    let base_coeff = 0.089 * a1 * 0.41 * lambda;

    // 3. 第一部分几何项: [ 2fi / (2fi - 1) ]^0.41
    let term1 = (2.0 * fi / (2.0 * fi - 1.0)).powf(0.41);

    // 4. 第二部分组合项: [ gamma^0.3 * (1 - gamma)^1.39 ] / [ (1 + gamma)^(1/3) ]
    let numerator_term2 = gamma.powf(0.3) * (1.0 - gamma).powf(1.39);
    let denominator_term2 = (1.0 + gamma).powf(1.0 / 3.0);
    let term2 = numerator_term2 / denominator_term2;

    // 5. 中括号内的嵌套复杂项（Dw 已完全消去）
    let inner_geom_ratio = (fi / fe) * ((2.0 * fe - 1.0) / (2.0 * fi - 1.0));
    let inner_power = inner_geom_ratio.powf(0.41) * ((1.0 - gamma) / (1.0 + gamma)).powf(1.72);

    // 6. 带有 10/3 次方和 -3/10 次方的外层括号项
    let bracket_base = 1.0 + (1.04 * inner_power).powf(10.0 / 3.0);
    let term3 = bracket_base.powf(-3.0 / 10.0);

    // 7. 最终组合
    Ok(base_coeff * term1 * term2 * term3)
}
// ISO 76 中静载荷系数
pub fn f0r(gamma: f64, fi: f64, fe: f64, sigma_max: f64) -> anyhow::Result<f64> {
    // 1. 计算f_rho
    let (_, f_rho_i) = calc_sum_rho_f_rho_i(gamma, 1.0, fi)?;
    // 1. 计算 e, kappa
    // 内圈
    let (_, _, _, e_i, kappa_i) = hertz_coefficients(f_rho_i)?;
    let m_i = e_i.powf(2.0);
    let term1 = sigma_max / 4000.0;
    let numerator_term2 = ellipe_agm(m_i)?;
    let denominator_term2 = 2.0 + gamma / (1.0 - gamma) - 1.0 / (2.0 * fi);
    let term2 = numerator_term2 / denominator_term2;
    let f0i = 2.072 * term1.powf(3.0) * kappa_i * term2.powf(2.0);

    // 外圈

    let (_, f_rho_e) = calc_sum_rho_f_rho_o(gamma, 1.0, fe)?;
    let (_, _, _, e_e, kappa_e) = hertz_coefficients(f_rho_e)?;
    let m_e = e_e.powf(2.0);
    let numerator_term2e = ellipe_agm(m_e)?;
    let denominator_term2e = 2.0 - gamma / (1.0 + gamma) - 1.0 / (2.0 * fe);
    let term2e = numerator_term2e / denominator_term2e;
    let f0e = 2.072 * term1.powf(3.0) * kappa_e * term2e.powf(2.0);
    // 取两者最小的值
    let f0 = f0e.min(f0i);
    Ok(f0)
}

#[cfg(test)]
mod tests {
    use anyhow::Ok;

    use crate::{
        assert_approx_eq,
        iso::{f0r, fcr},
    };

    #[test]
    fn test_fcr_iso281() -> anyhow::Result<()> {
        // 测试数据来自 GB/T 6391-2010 表5.1.2.2
        let eps = 1e-1;
        let a1 = 98.0665 / 0.089;
        let lambda = 0.95;
        let fi = 0.52;
        let fe = 0.52;
        let gamma = 0.4;
        let fc = fcr(a1, lambda, fi, fe, gamma)?;
        assert_approx_eq!(fc, 48.4, eps);
        let gamma = 0.3;
        let fc = fcr(a1, lambda, fi, fe, gamma)?;
        assert_approx_eq!(fc, 56.0, eps);
        let gamma = 0.2;
        let fc = fcr(a1, lambda, fi, fe, gamma)?;
        assert_approx_eq!(fc, 59.9, eps);
        let gamma = 0.1;
        let fc = fcr(a1, lambda, fi, fe, gamma)?;
        assert_approx_eq!(fc, 55.5, eps);
        let gamma = 0.01;
        let fc = fcr(a1, lambda, fi, fe, gamma)?;
        assert_approx_eq!(fc, 29.1, eps);
        Ok(())
    }
    #[test]
    #[ignore]
    fn test_print_fcr_iso281() -> anyhow::Result<()> {
        let a1 = 98.0665 / 0.089;
        let lambda = 0.95;
        let fi = 0.52;
        let fe = 0.52;
        let mut gamma = 0.01;
        println!("{}\t{}", "gamma", "fc");
        while gamma < 0.4 + 1e-5 {
            let fc = fcr(a1, lambda, fi, fe, gamma)?;
            println!("{:.2}\t{:.1}", gamma, fc);
            gamma = gamma + 0.01;
        }
        Ok(())
    }
    #[test]
    fn test_f0r_iso76() -> anyhow::Result<()> {
        let eps = 0.1;
        let fi = 0.52;
        let fe = 0.53;
        let sigma_max = 4200.0;
        let gamma = 0.1;
        let f0 = f0r(gamma, fi, fe, sigma_max)?;
        assert_approx_eq!(f0, 16.4, eps);
        let gamma = 0.2;
        let f0 = f0r(gamma, fi, fe, sigma_max)?;
        assert_approx_eq!(f0, 14.0, eps);
        let gamma = 0.3;
        let f0 = f0r(gamma, fi, fe, sigma_max)?;
        assert_approx_eq!(f0, 11.6, eps);
        let gamma = 0.4;
        let f0 = f0r(gamma, fi, fe, sigma_max)?;
        assert_approx_eq!(f0, 9.4, eps);
        let gamma = 0.01;
        let f0 = f0r(gamma, fi, fe, sigma_max)?;
        assert_approx_eq!(f0, 14.7, eps);
        Ok(())
    }
    #[test]
    #[ignore]
    fn test_print_f0r_iso76() -> anyhow::Result<()> {
        let fi = 0.52;
        let fe = 0.53;
        let sigma_max = 4200.0;
        let mut gamma = 0.00;
        println!("{}\t{}", "gamma", "f0");
        while gamma < 0.4 + 1e-5 {
            let f0 = f0r(gamma, fi, fe, sigma_max)?;
            println!("{:.2}\t{:.1}", gamma, f0);
            gamma += 0.01;
        }
        Ok(())
    }
}
