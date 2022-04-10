use thiserror::Error;
use std::path::PathBuf;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },

    #[error("failed solving the problem")]
    SolverProblem,
}