use std::f64::consts::PI;

use roots::{SimpleConvergency, find_root_brent};

use crate::common::{hertz_coefficients, jra};

// 计算输入
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
pub struct ACBOut {
    pub cr: f64,
    pub c0r: f64,
    pub ca: f64,
    pub c0a: f64,
}

// 载荷输出结果
pub struct LoadResult {
    pub psi: Vec<f64>,
    pub load: Vec<f64>,
    pub sigma: Vec<f64>,
}
// 计算载荷分布函数
// 返回一个关于方位角psi的函数
pub fn load_distribution(
    fr: f64,
    fa: f64,
    alpha: f64,
    z: i32,
) -> anyhow::Result<impl Fn(f64) -> f64> {
    let helper = |epsilon| -> f64 {
        let (jr, ja) = jra(epsilon, 1.5).unwrap();
        fr * alpha.tan() / fa - jr / ja
    };
    let mut convergency = SimpleConvergency {
        eps: 1e-15,
        max_iter: 100,
    };
    let epsilon = find_root_brent(1e-5, 1e5 - 1e-5, helper, &mut convergency)?;
    // println!("{}", epsilon);
    let (jr, _ja) = jra(epsilon, 1.5)?;
    let q_max = fr / (jr * z as f64 * alpha.sin());
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

// 计算应力
pub fn calc_sigma(q: f64, a_star: f64, b_star: f64, sigma_rho: f64) -> anyhow::Result<f64> {
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
// 计算载荷分布-完整代码
// pub fn calc_dist_all(acb_input: ACBInput) -> anyhow::Result<LoadResult> {
//     let z = acb_input.z;
//     let mut psi_vec = vec![0f64; z as usize];
//     let mut load_vec = vec![0f64; z as usize];
//     let mut sigma_vec = vec![0f64; z as usize];
//     let load_fun = load_distribution(acb_input.fr, acb_input.fa, acb_input.alpha, acb_input.z)?;
//     let (sigma_rho, f_rho) = ca();
//     let (a_star, b_star, _, _, _) = hertz_coefficients(f_rho)?;
//     for i in 0..z as usize {
//         psi_vec[i] = i as f64 * 2.0 * PI / z as f64;
//         load_vec[i] = load_fun(psi_vec[i]);
//         sigma_vec[i] = calc_sigma(load_vec[i], a_star, b_star, sigma_rho)?;
//     }
//     let load_result = LoadResult {
//         psi: psi_vec,
//         load: load_vec,
//         sigma: sigma_vec,
//     };

//     Ok(load_result)
// }
#[cfg(test)]
mod tests {

    use super::*;
    use crate::macros::assert_approx_eq;
    #[test]
    fn test() -> anyhow::Result<()> {
        let load_fun = load_distribution(0.6, 1.0, PI / 4.0, 16)?;
        let res = load_fun(0.0);
        println!("{}", res);
        Ok(())
    }
    
    #[test]
    fn test_rho() -> anyhow::Result<()> {
        // 根据《滚动轴承分析》第一卷例2.5的测试用例
        let gamma = 0.1954;
        let fi = 0.52;
        let fo = 0.52;
        let d = 12.7;
        let (sum_rho_i, f_rho_i) = calc_sum_rho_f_rho_i(gamma, d, fi)?;
        let (sum_rho_o, f_rho_o) = calc_sum_rho_f_rho_o(gamma, d, fo)?;
        assert_approx_eq!(sum_rho_i, 0.202);
        assert_approx_eq!(f_rho_i, 0.9399);
        assert_approx_eq!(sum_rho_o, 0.1378);
        assert_approx_eq!(f_rho_o, 0.9120);
        Ok(())
    }
}
