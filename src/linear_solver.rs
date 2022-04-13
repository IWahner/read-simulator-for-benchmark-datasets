extern crate lp_modeler;

use lp_modeler::solvers::{CbcSolver, SolverTrait};
use lp_modeler::dsl::*;
use lp_modeler::constraint;
use lp_modeler::dsl::variables::lp_sum;
use std::collections::HashMap;
use rand::seq::SliceRandom;
use rand::rngs::StdRng;
use rand::SeedableRng;


pub struct LinearSolver{
}

impl Default for LinearSolver {
    fn default() -> Self {
         Self::new()
      }
}
 
impl LinearSolver{
 
     pub fn new() -> Self {
         LinearSolver {
         }
    }

    pub fn variantselection(&self, genloci: i32, mut numberofvariants: i32, wanted_freq: f32, seed: u64) -> HashMap<String, f32> {
        let mut sol = HashMap::new();
        let mut beginning: i32 = 0;
        let for_calculation = numberofvariants as f32; 
            while numberofvariants > 0 {   
            //tested, solvers works with 20000 variables
            let partialnumber: i32 = if numberofvariants >= 20000 {
                20000
            }
            else {
                numberofvariants
            };
            //let fraction = partialnumber as f32;
            let partial_genloci = genloci as f32 *(partialnumber as f32/for_calculation);
            println!("partial_gen: {}", partial_genloci);
            numberofvariants -= 20000;
            let mut v = Vec::with_capacity(partialnumber.try_into().unwrap());
            let end = beginning + partialnumber;
            //Vector of binaries
            for i in beginning..end {
                v.push(LpBinary::new(&("x".to_owned()+&i.to_string())));
            }
            let mut rng = StdRng::seed_from_u64(seed);
            v.shuffle(&mut rng);
    
            beginning += partialnumber;
            //println!("{:?}", v); //Inhalt des Vectors
    
            let ref actual_freq = LpContinuous::new("d");
    
            // Define problem and objective sense
            let mut problem = LpProblem::new("One Problem", LpObjective::Minimize);
    
            // Objective Function: Minimize ds - Ds
            problem += actual_freq - wanted_freq;
    
            // Constraint:
            problem += lp_sum(&v).equal(actual_freq*partial_genloci);
            
            problem += constraint!(actual_freq - wanted_freq >= 0);
        
            // Specify solver
            let solver = CbcSolver::new();
    
            // Run optimisation and process output vcf and bam file
            match solver.run(&problem) {
                Ok(solution) => {
                    println!("Status {:?}", solution.status);
                    sol.extend(solution.results)
                },
                Err(msg) => println!("{}", msg),
            }
        }
        sol
    }
}