use pyo3_stub_gen::Result;

fn main() -> Result<()> {
    // `stub_info` is a function defined by `define_stub_info_gatherer!` macro.
    // $env:PATH = "$(python -c 'import sys; print(sys.base_prefix)');$env:PATH"
    let stub = python_wrapper::stub_info()?;
    stub.generate()?;
    Ok(())
}
