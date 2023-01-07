use pyo3::prelude::*;
use std::path::PathBuf;
use relate::pipelines::make_chunks as make_chunks_rs;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn make_chunks(haps: PathBuf, sample: PathBuf, map: PathBuf, output: PathBuf, dist: Option<PathBuf>, use_transitions: bool, memory: f32) -> PyResult<()> {
    Ok(make_chunks_rs(haps, sample, map, dist.unwrap_or(PathBuf::from("unspecified")), output, !use_transitions, memory).unwrap())
}

/// A Python module implemented in Rust.
#[pymodule]
fn relatepy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(make_chunks, m)?)?;
    Ok(())
}