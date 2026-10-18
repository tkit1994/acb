use rustmex::prelude::*;

// NOTE: these might move into the prelude some day
use rustmex::{
    // Trait describing 'Matlab Classes', such as the aforementioned Numeric array.
    convert::{IntoRust, ToMatlab},
};

use acb_core::iso::f0;

#[rustmex::entrypoint]
fn calc_f0(lhs: Lhs, rhs: Rhs) -> rustmex::Result<()> {
    let coefficient = rhs.get(0).unwrap().to_rust()?;
    let gamma = rhs.get(1).unwrap().to_rust()?;
    let fi = rhs.get(2).unwrap().to_rust()?;
    let fe = rhs.get(3).unwrap().to_rust()?;
    let sigma_max = rhs.get(4).unwrap().to_rust()?;
    let res = f0(coefficient, gamma, fi, fe, sigma_max).unwrap();
    let ret = lhs.get_mut(0).unwrap();
    ret.replace(res.to_matlab());
    Ok(())
}
