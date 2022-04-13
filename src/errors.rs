use thiserror::Error;
use std::path::PathBuf;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    #[error("file not found: {path}")]
    FileNotFound { path: PathBuf },

    #[error("failed solving the problem")]
    SolvingProblem,

    #[error("wrong input for solver")]
    SolverInputProblem,
}