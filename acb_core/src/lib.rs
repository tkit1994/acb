use anyhow::Ok;

use crate::core::ACBInput;
use crate::core::ACBOut;
use crate::core::calc_c0a;
use crate::core::calc_c0r;
use crate::core::calc_ca;
use crate::core::calc_cr;
use crate::core::calc_dist_all;
use crate::core::calc_life_lp;
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
    let l10 = calc_life_lp(acb_input, None)?;
    let acb_out = ACBOut {
        cr,
        c0r,
        ca,
        c0a,
        load_result,
        p0a,
        p0r,
        s0,
        l10,
    };
    Ok(acb_out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn test_calc_all() -> anyhow::Result<()> {
        let acb_input = ACBInput::default();
        let acb_out = calc_all(acb_input)?;
        println!("{:#.4?}", acb_out);

        Ok(())
    }
    #[test]
    #[ignore]
    fn test_calc_all_1() -> anyhow::Result<()> {
        let acb_input = ACBInput{
            dw: 7.144,
            dpw: 45.0,
            z: 18,
            alpha: 45f64.to_radians(),
            fr: 1.0,
            fa: 5000.0,
            fi: 0.54,
            fe: 0.54,
            bm: 1.3,
            lambda: 0.95,
        };
        let acb_out = calc_all(acb_input)?;
        println!("{:#.4?}", acb_out);

        Ok(())
    }
}
