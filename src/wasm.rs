use serde_wasm_bindgen::{from_value, to_value};
use wasm_bindgen::prelude::*;

use crate::{calc_all, core::ACBInput};

#[wasm_bindgen]
pub fn calc(json: JsValue) -> Result<JsValue, JsValue> {
    let acb_input: ACBInput = from_value(json)?;
    let res = calc_all(acb_input).map_err(|err| JsValue::from_str(&format!("{}", err)))?;
    let out_json = to_value(&res)?;
    Ok(out_json)
}
