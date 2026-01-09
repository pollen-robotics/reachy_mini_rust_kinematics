#[cfg(feature = "python")]
use pyo3_stub_gen::Result;
#[cfg(feature = "python")]
use reachy_mini_rust_kinematics::python::stub_info;

#[cfg(feature = "python")]
fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().filter_or("RUST_LOG", "info")).init();
    let stub = stub_info()?;
    stub.generate()?;
    Ok(())
}
