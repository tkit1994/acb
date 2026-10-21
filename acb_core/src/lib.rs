use anyhow::Ok;

use crate::core::ACBInput;
use crate::core::ACBOut;
use crate::core::calc_c0a;
use crate::core::calc_c0r;
use crate::core::calc_ca;
use crate::core::calc_cr;
use crate::core::calc_dist_all;
use crate::core::calc_p0a;
use crate::core::calc_p0r;

pub mod common;
pub mod core;
pub mod iso;
pub mod macros;

pub fn calc_all(acb_input: ACBInput) -> anyhow::Result<ACBOut> {
    let c0a = calc_c0a(acb_input)?;
    let c0r = calc_c0r(acb_input)?;
    let ca = calc_ca(acb_input)?;
    let cr = calc_cr(acb_input)?;
    let load_result = calc_dist_all(acb_input)?;
    let p0r = calc_p0r(acb_input)?;
    let p0a = calc_p0a(acb_input)?;
    let s0 = (c0r / p0r).min(c0a / p0a);
    let acb_out = ACBOut {
        cr,
        c0r,
        ca,
        c0a,
        load_result,
        p0a,
        p0r,
        s0,
    };
    Ok(acb_out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn test_calc_all() -> anyhow::Result<()> {
        let acb_input = ACBInput {
            dw: 22.23,
            dpw: 125.3,
            z: 16,
            alpha: 40.0f64.to_radians(),
            fr: 17800.0,
            fa: 17800.0,
            fi: 11.6 / 22.23,
            fe: 11.6 / 22.23,
            bm: 1.0,
        };
        let acb_out = calc_all(acb_input)?;
        println!("{:#.4?}", acb_out);

        Ok(())
    }
}
