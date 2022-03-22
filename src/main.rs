extern crate lp_modeler;
//use crate::rust_htslib::bcf::{Reader, Read};

use lp_modeler::solvers::{CbcSolver, SolverTrait};
use lp_modeler::dsl::*;
use lp_modeler::constraint;
use lp_modeler::dsl::variables::lp_sum;
use lp_modeler::solvers::Solution;
use rust_htslib::bcf::Reader;
use rust_htslib::bcf::Read;
use std::convert::TryFrom;


fn main() {
    // Define problem variables

    //Number of Loci from the vcf containing all Entries of Chr21
    let path = &"Rohdaten/only21.vcf";
    let mut all21 = Reader::from_path(path).expect("Error opening file.");
    let genloci = all21.records().count() as i32;
    println!("{}", genloci);

    //Numberofvariants from the vcf filtered by the bed file
    /*let second_path = &"Rohdaten/filtered.vcf";
    let mut filtered = Reader::from_path(second_path).expect("Error opening file.");
    let variants_count = filtered.records().count();
    let mut v = Vec::with_capacity(variants_count);
    let numberofvariants = variants_count as i32;
    println! ("{}", numberofvariants);*/
    //let numberofvariants = 57648;
    let numberofvariants = 5000;    
    let mut v = Vec::with_capacity(numberofvariants);
    let wanted_freq = 0.03; // gewuenschte Variantenfrequenz

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
    
    problem += constraint!(actual_freq - wanted_freq >= 0);
 
    // Specify solver
    let solver = CbcSolver::new();

    // Run optimisation and process output hashmap
    match solver.run(&problem) {
        Ok(solution) => {
            variantselection(solution);
            /*println!("Status {:?}", solution.status);
            for (name, value) in solution.results.iter() {
                println!("value of {} = {}", name, value);
            }*/
        },
        Err(msg) => println!("{}", msg),
    }
}

fn variantselection(solution: Solution){
    println!("in der Funktion");
    println!("Status {:?}", solution.status);
    for (name, value) in solution.results.iter() {
        //println!("value of {} = {}", name, value);
    } 
}
