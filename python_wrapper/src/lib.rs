use pyo3::prelude::*;
use pyo3_stub_gen::define_stub_info_gatherer;

/// A Python module implemented in Rust.
#[pymodule]
mod python_wrapper {
    use acb_core::calc_all;
    use acb_core::core::{ACBInput, ACBOut};
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};
    #[gen_stub_pyclass]
    #[pyclass(name = "ACBInput")]
    struct PyACBInput {
        pub inner: ACBInput,
    }
    #[gen_stub_pymethods]
    #[pymethods]
    impl PyACBInput {
        #[new]
        fn new() -> Self {
            let acb_input = ACBInput::default();
            Self { inner: acb_input }
        }
    }
    #[gen_stub_pyclass]
    #[pyclass(name = "ACBOut")]
    struct PyACBOut {
        pub inner: ACBOut,
    }
    #[gen_stub_pymethods]
    #[pymethods]
    impl PyACBOut {
        fn __repr__(&self) -> String {
            format!("ACBOut: {:?}", self.inner)
        }
        fn getload(&self) -> Vec<f64> {
            self.inner.load_result.load.to_owned()
        }
    }
    #[gen_stub_pyfunction]
    #[pyfunction]
    fn calc(input: &PyACBInput) -> PyResult<PyACBOut> {
        let inner = calc_all(input.inner).map_err(|err| PyValueError::new_err(err.to_string()))?;
        let res = PyACBOut { inner };
        Ok(res)
    }
}

define_stub_info_gatherer!(stub_info);
