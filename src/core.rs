use std::f64::consts::PI;

use roots::{SimpleConvergency, find_root_brent};

use crate::common::{hertz_coefficients, jra};

// 计算输入
#[derive(Debug)]
pub struct ACBInput {
    pub dw: f64,
    pub dpw: f64,
    pub z: i32,
    pub alpha: f64,
    pub fr: f64,
    pub fa: f64,
    pub fi: f64,
    pub fe: f64,
}

// 计算输出结果
#[derive(Debug)]
pub struct ACBOut {
    pub cr: f64,
    pub c0r: f64,
    pub ca: f64,
    pub c0a: f64,
}

// 载荷输出结果
#[derive(Debug)]
pub struct LoadResult {
    pub psi: Vec<f64>,
    pub load: Vec<f64>,
    pub sigma_i: Vec<f64>,
    pub sigma_o: Vec<f64>,
}

// 1. 根据 fr, fa, alpha 计算 epsilon
pub fn calc_epsilon(fr: f64, fa: f64, alpha: f64) -> anyhow::Result<f64> {
    let helper = |epsilon| -> f64 {
        let (jr, ja) = jra(epsilon, 1.5).unwrap();
        fr * alpha.tan() / fa - jr / ja
    };
    let mut convergency = SimpleConvergency {
        eps: 1e-15,
        max_iter: 100,
    };
    let epsilon = find_root_brent(1e-5, 1e5 - 1e-5, helper, &mut convergency)?;
    Ok(epsilon)
}
// 2. 计算q_max
pub fn calc_q_max(epsilon: f64, z: f64, alpha: f64, fr: f64) -> anyhow::Result<f64> {
    let (jr, _ja) = jra(epsilon, 1.5)?;
    let q_max = fr / (jr * z * alpha.cos());
    Ok(q_max)
}
// 3. 计算载荷分布函数
// 返回一个关于方位角psi的函数, 函数返回值是方位角psi出钢球受力
pub fn load_distribution(
    epsilon: f64,
    z: f64,
    alpha: f64,
    fr: f64,
) -> anyhow::Result<impl Fn(f64) -> f64> {
    let q_max = calc_q_max(epsilon, z, alpha, fr)?;
    let load_fun = move |psi: f64| -> f64 {
        let term = 1.0 - 1.0 / (2.0 * epsilon) * (1.0 - psi.cos());
        if term < 0.0 {
            0.0
        } else {
            q_max * term.powf(1.5)
        }
    };
    Ok(load_fun)
}
// 4.1
pub fn calc_sigma(q: f64, a_star: f64, b_star: f64, sigma_rho: f64) -> anyhow::Result<f64> {
    if q.abs() < 1e-3 {
        return Ok(0.0);
    } 
    let term = (q / sigma_rho).powf(1.0 / 3.0);
    let a = 0.0236 * a_star * term;
    let b = 0.0236 * b_star * term;
    let sigma = 3.0 * q / (2.0 * PI * a * b);
    Ok(sigma)
}

// 计算曲率函数, 曲率和
// sum_rho, f_rho
pub fn calc_sum_rho_f_rho_i(gamma: f64, d: f64, fi: f64) -> anyhow::Result<(f64, f64)> {
    let sum_rho = 1.0 / d * (4.0 - 1.0 / fi + 2.0 * gamma / (1.0 - gamma));
    let f_rho =
        (1.0 / fi + 2.0 * gamma / (1.0 - gamma)) / (4.0 - 1.0 / fi + 2.0 * gamma / (1.0 - gamma));
    Ok((sum_rho, f_rho))
}

pub fn calc_sum_rho_f_rho_o(gamma: f64, d: f64, fo: f64) -> anyhow::Result<(f64, f64)> {
    let sum_rho = 1.0 / d * (4.0 - 1.0 / fo - 2.0 * gamma / (1.0 + gamma));
    let f_rho =
        (1.0 / fo - 2.0 * gamma / (1.0 + gamma)) / (4.0 - 1.0 / fo - 2.0 * gamma / (1.0 + gamma));
    Ok((sum_rho, f_rho))
}

pub fn deg2rad(alpha: f64) -> f64 {
    alpha * PI / 180.0
}

pub fn rad2deg(alpha: f64) -> f64 {
    alpha * 180.0 / PI
}
// 计算载荷分布-完整代码
pub fn calc_dist_all(acb_input: ACBInput) -> anyhow::Result<LoadResult> {
    let fr = acb_input.fr;
    let fa = acb_input.fa;
    let alpha = acb_input.alpha;
    let z = acb_input.z as f64;
    let dw = acb_input.dw;
    let dm = acb_input.dpw;
    let fi = acb_input.fi;
    let fo = acb_input.fe;
    // 1. 计算epsilon
    let epsilon = calc_epsilon(fr, fa, alpha)?;
    // 2. 计算载荷分布
    let load_fun = load_distribution(epsilon, z, alpha, fr)?;
    // 3. 计算应力
    // 3.1 计算曲率和，曲率差, hertz
    let gamma = dw * alpha.cos() / dm;
    // 内圈
    let (sum_rho_i, f_rho_i) = calc_sum_rho_f_rho_i(gamma, dw, fi)?;
    // 外圈
    let (sum_rho_o, f_rho_o) = calc_sum_rho_f_rho_o(gamma, dw, fo)?;
    // 3.2 Hertz 接触系数
    // 内圈
    let (a_star_i, b_star_i, _, _, _) = hertz_coefficients(f_rho_i)?;
    // 外圈
    let (a_star_o, b_star_o, _, _, _) = hertz_coefficients(f_rho_o)?;

    // 3.3 计算载荷分布
    let mut psi_vec = vec![0f64; z as usize];
    let mut load_vec = vec![0f64; z as usize];
    let mut sigma_vec_i = vec![0f64; z as usize];
    let mut sigma_vec_o = vec![0f64; z as usize];
    for i in 0..z as usize {
        psi_vec[i] = i as f64 * 2.0 * PI / z;
        load_vec[i] = load_fun(psi_vec[i]);
        sigma_vec_i[i] = calc_sigma(load_vec[i], a_star_i, b_star_i, sum_rho_i)?;
        sigma_vec_o[i] = calc_sigma(load_vec[i], a_star_o, b_star_o, sum_rho_o)?;
    }
    let load_result = LoadResult {
        psi: psi_vec,
        load: load_vec,
        sigma_i: sigma_vec_i,
        sigma_o: sigma_vec_o,
    };

    Ok(load_result)
}
#[cfg(test)]
mod tests {

    use anyhow::Ok;

    use super::*;
    use crate::assert_approx_eq;
    #[test]
    fn test_rho() -> anyhow::Result<()> {
        // 根据《滚动轴承分析》第一卷例2.5的测试用例
        let gamma = 0.1954;
        let fi = 0.52;
        let fo = 0.52;
        let dw = 12.7;
        let (sum_rho_i, f_rho_i) = calc_sum_rho_f_rho_i(gamma, dw, fi)?;
        let (sum_rho_o, f_rho_o) = calc_sum_rho_f_rho_o(gamma, dw, fo)?;
        assert_approx_eq!(sum_rho_i, 0.202);
        assert_approx_eq!(f_rho_i, 0.9399);
        assert_approx_eq!(sum_rho_o, 0.1378);
        assert_approx_eq!(f_rho_o, 0.9120);
        Ok(())
    }

    #[test]
    fn test_epsilon() -> anyhow::Result<()> {
        // 测试数据来源 滚动轴承分析 例7.7
        let fr = 17800.0;
        let fa = 17800.0;
        let alpha = deg2rad(40.0);
        let epsilon = calc_epsilon(fr, fa, alpha)?;
        assert_approx_eq!(epsilon, 0.455, 2e-3);
        Ok(())
    }

    #[test]
    fn test_q_max() -> anyhow::Result<()> {
        // 测试数据来源 滚动轴承分析 例7.7
        let epsilon = 0.455;
        let z = 16.0;
        let alpha = deg2rad(40.0);
        let fr = 17800.0;
        let q_max = calc_q_max(epsilon, z, alpha, fr)?;
        assert_approx_eq!(q_max, 6571.0, 2.0);
        Ok(())
    }
    #[test]
    fn test_load_distribution() -> anyhow::Result<()> {
        // 测试数据来源 滚动轴承分析 例7.7
        let eps = 10.0;
        let epsilon = 0.455;
        let z = 16.0;
        let alpha = deg2rad(40.0);
        let fr = 17800.0;
        let load_fun = load_distribution(epsilon, z, alpha, fr)?;
        assert_approx_eq!(6571.0, load_fun(deg2rad(0.0)), eps);
        assert_approx_eq!(5765.0, load_fun(deg2rad(22.5)), eps);
        assert_approx_eq!(3670.0, load_fun(deg2rad(45.0)), eps);
        assert_approx_eq!(1200.0, load_fun(deg2rad(67.5)), eps);
        assert_approx_eq!(0.0, load_fun(deg2rad(90.0)), eps);
        Ok(())
    }
    #[test]
    fn test_sigma() -> anyhow::Result<()> {
        // 测试数据来源 滚动轴承分析 例6.1
        let dw = 22.23;
        let dm = 125.3;
        let alpha = deg2rad(38.9);
        let fi = 11.63 / 22.23;

        // gamma
        let gamma = dw * alpha.cos() / dm;
        let (sum_rho, f_rho) = calc_sum_rho_f_rho_i(gamma, dw, fi)?;
        let (a_star, b_star, _, _, _) = hertz_coefficients(f_rho)?;
        let q = 3543.0;
        let sigma = calc_sigma(q, a_star, b_star, sum_rho)?;
        assert_approx_eq!(1976.0, sigma, 10.0);
        Ok(())
    }
    #[test]
    #[ignore]
    fn test_load_distribution_all() -> anyhow::Result<()> {
        let acb_input = ACBInput {
            dw: 22.23,
            dpw: 125.3,
            z: 16,
            alpha: deg2rad(40.0),
            fr: 17800.0,
            fa: 17800.0,
            fi: 11.6 / 22.23,
            fe: 11.6 / 22.23,
        };
        let res = calc_dist_all(acb_input)?;
        println!("{:#.4?}", res);
        Ok(())
    }
}
