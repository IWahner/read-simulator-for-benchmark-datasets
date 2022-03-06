extern crate lp_modeler;

use lp_modeler::solvers::{CbcSolver, SolverTrait};
use lp_modeler::dsl::*;
//use lp_modeler::constraint;
use lp_modeler::dsl::variables::lp_sum;


fn main() {
    // Define problem variables
    let wanted_freq = 0.3; // gewuenschte Variantenfrequenz
    let genloci = 21; //Anzahl der Genloci
    let numberofvariants = 15;
    let mut v = Vec::with_capacity(numberofvariants);

    //Vector der Binäries
    for i in 0..numberofvariants {
        v.push(LpBinary::new(&("x".to_owned()+&i.to_string())));
    }

    //println!("{:?}", v); //Inhalt des Vectors

    let ref actual_freq = LpContinuous::new("d");
    // Define problem and objective sense
    let mut problem = LpProblem::new("One Problem", LpObjective::Minimize);

    // Objective Function: Minimize ds - Ds
    problem += actual_freq - wanted_freq;

    // Constraint:
    problem += lp_sum(&v).equal(actual_freq*genloci);

    // Specify solver
    let solver = CbcSolver::new();

    // Run optimisation and process output hashmap
    match solver.run(&problem) {
        Ok(solution) => {
            println!("Status {:?}", solution.status);
            for (name, value) in solution.results.iter() {
                println!("value of {} = {}", name, value);
            }
        },
        Err(msg) => println!("{}", msg),
    }
}
