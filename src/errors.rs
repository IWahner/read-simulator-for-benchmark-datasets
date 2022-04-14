use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum Error {
    #[error("file not found")]
    FileNotFound,

    #[error("failed solving the problem")]
    SolvingProblem,

    #[error("wrong input for solver")]
    SolverInputProblem,

    #[error("invalid input for wanted frequency")]
    InvalidFrequency,
}