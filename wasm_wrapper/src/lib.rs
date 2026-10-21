use std::fmt::Display;

use serde_wasm_bindgen::{from_value, to_value};
use wasm_bindgen::prelude::*;

use acb_core::{
    calc_all,
    core::{ACBInput, calc_epsilon, load_distribution},
};

fn err2js(err: impl Display) -> JsValue {
    JsValue::from_str(&format!("{}", err))
}
#[wasm_bindgen]
pub struct LoadDistribution {
    load_fn: Box<dyn Fn(f64) -> f64>,
}
#[wasm_bindgen]
impl LoadDistribution {
    #[wasm_bindgen(constructor)]
    pub fn new(json: JsValue) -> Result<LoadDistribution, JsValue> {
        let acb_input: ACBInput = from_value(json)?;
        let z = acb_input.z as f64;
        let alpha = acb_input.alpha;
        let fr = acb_input.fr;
        let fa = acb_input.fa;

        let epsilon = calc_epsilon(fr, fa, alpha).map_err(err2js)?;
        let load_fn = load_distribution(epsilon, z, alpha, fr).map_err(err2js)?;
        Ok(Self {
            load_fn: Box::new(load_fn),
        })
    }
    #[wasm_bindgen]
    pub fn load(&self, psi: f64) -> Result<f64, JsValue> {
        let res = (self.load_fn)(psi);
        Ok(res)
    }
}

#[wasm_bindgen]
pub fn calc(json: JsValue) -> Result<JsValue, JsValue> {
    let acb_input: ACBInput = from_value(json)?;
    let res = calc_all(acb_input).map_err(|err| JsValue::from_str(&format!("{}", err)))?;
    let out_json = to_value(&res)?;
    Ok(out_json)
}
