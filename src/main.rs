extern crate lp_modeler;

use lp_modeler::solvers::{CbcSolver, SolverTrait};
use lp_modeler::dsl::*;
use lp_modeler::constraint;
use lp_modeler::dsl::variables::lp_sum;
use lp_modeler::solvers::Solution;
use rust_htslib::bcf::*;
use rand::seq::SliceRandom;
use rand::rngs::StdRng;
use rand::SeedableRng;
use code::output_writer::OutputWriter;
use std::time::{Instant};


fn main() {
    // Define problem variables
    let now = Instant::now();
    //Number of Loci equals the full length of chr21
    let path = &"Rohdaten/filtered.vcf";
    let all21 = Reader::from_path(path).expect("Error opening file.");
    let allheader = all21.header().header_records();
    let mut genloci = 0;
    for hdr_records in allheader {
       genloci = match hdr_records {
            HeaderRecord::Contig{values, ..} => values.get("length").unwrap().parse().unwrap(),
            _ => continue,
        };

    }
    
    /*//Numberofvariants from the vcf for CHM1
    let second_path = &"Rohdaten/only21.vcf";
    let mut numberofvariants = 0;
    let mut to_filter = Reader::from_path(second_path).expect("Error opening file.");
    for gta in to_filter.records(){
        let record = gta.expect("Fail to read record");
        let gts = record.genotypes().expect("Error reading genotypes");
        let geno1 = gts.get(0);
        match geno1[0].index(){
            Some(1) => numberofvariants += 1,
            Some(_) => continue,
            None => continue
        };

    }*/
    println! ("Number of genloci {}", genloci);   
    //let numberofvariants = genloci;  //Stackoverflow with this Number of variants
    let numberofvariants = 10000;
    let mut v = Vec::with_capacity(numberofvariants);
    let wanted_freq = 0.0001; // gewuenschte Variantenfrequenz

    //Vector der Binäries
    for i in 0..numberofvariants {
        v.push(LpBinary::new(&("x".to_owned()+&i.to_string())));
    }
    let mut rng = StdRng::seed_from_u64(1337);
    v.shuffle(&mut rng);

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

    // Run optimisation and process output vcf and bam file
    match solver.run(&problem) {
        Ok(solution) => {
            //variantselection(solution);
            OutputWriter::new().variantselection(solution);
        },
        Err(msg) => println!("{}", msg),
    }
    let average_time = now.elapsed().as_millis();
    println!("average_time {} ms", average_time);
}






