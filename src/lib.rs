use pyo3::prelude::*;
use std::path::PathBuf;
use relate::pipelines::{MakeChunks, Paint};

/// Formats the sum of two numbers as string.
#[pyfunction]
fn make_chunks(haps: PathBuf, sample: PathBuf, map: PathBuf, output: PathBuf, dist: Option<PathBuf>, use_transitions: Option<bool>, memory: Option<f32>) -> PyResult<()> {
    let options = MakeChunks::new(
        haps, sample, map, dist.unwrap_or(PathBuf::from("unspecified")), output, !use_transitions.unwrap_or(true), memory.unwrap_or(5.)
    );
    Ok(options.execute().unwrap())
}

#[pyfunction]
fn paint(output: PathBuf, chunk_index: usize, theta: f64, rho: f64) -> PyResult<()> {
    let painting = vec![theta, rho];
    let options = Paint::new(chunk_index, output, painting);
    Ok(options.execute().unwrap())
}

/// A Python module implemented in Rust.
#[pymodule]
fn relatepy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(make_chunks, m)?)?;
    m.add_function(wrap_pyfunction!(paint, m)?)?;
    Ok(())
}