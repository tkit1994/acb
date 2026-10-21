use pyo3::prelude::*;
/// A Python module implemented in Rust.
#[pymodule]
mod python_wrapper {
    use acb_core::calc_all;
    use acb_core::core::{ACBInput, ACBOut};
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;

    #[pyclass(name = "ACBInput")]
    struct PyACBInput {
        pub inner: ACBInput,
    }

    #[pymethods]
    impl PyACBInput {
        #[new]
        fn new() -> Self {
            let acb_input = ACBInput {
                dw: 22.23,
                dpw: 125.3,
                z: 16,
                alpha: 40.0f64.to_radians(),
                fr: 17800.0,
                fa: 17800.0,
                fi: 11.6 / 22.23,
                fe: 11.6 / 22.23,
            };
            Self { inner: acb_input }
        }
    }
    #[pyclass(name = "ACBOut")]
    struct PyACBOut {
        pub inner: ACBOut,
    }
    #[pymethods]
    impl PyACBOut {
        fn __repr__(&self) -> String {
            format!("ACBOut: {:?}", self.inner)
        }
        fn getload(&self) -> Vec<f64> {
            self.inner.load_result.load.to_owned()
        }
    }
    #[pyfunction]
    fn calc(input: &PyACBInput) -> PyResult<PyACBOut> {
        let inner = calc_all(input.inner).map_err(|err| PyValueError::new_err(err.to_string()))?;
        let res = PyACBOut { inner };
        Ok(res)
    }
}
