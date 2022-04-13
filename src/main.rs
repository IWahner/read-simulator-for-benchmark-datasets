extern crate lp_modeler;

use lp_modeler::solvers::{CbcSolver, SolverTrait};
use lp_modeler::dsl::*;
use lp_modeler::constraint;
use lp_modeler::dsl::variables::lp_sum;
use rust_htslib::bcf::*;
use rand::seq::SliceRandom;
use rand::rngs::StdRng;
use rand::SeedableRng;
use code::output_writer::OutputWriter;
//use code::errors::Error;
use code::greedy_solver::GreedySolver;
use std::time::{Instant};
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::Read as OtherRead;
use std::collections::HashMap;


#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct SimData {
    variantfreq: f32,
    rng_seed: u64,
}


fn main() {
    //time measurment
    let now = Instant::now();

    //reading from the input yaml
    let filename = &"Rohdaten/test.yml";
    let mut wanted_freq: f32 = 0.0;
    let mut seed: u64 = 0;
    match File::open(filename) {
        Ok(mut file) => {
            let mut content = String::new();
            file.read_to_string(&mut content).unwrap();
            let input_yml: SimData = serde_yaml::from_str(&content).unwrap();
            wanted_freq = input_yml.variantfreq;
            seed = input_yml.rng_seed;
        }
        Err(error) => {
            println!("There is an error {}: {}", filename, error);
        }
    }
    println!("{}", wanted_freq);
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
      
    //Numberofvariants from the vcf for CHM1
    let second_path = &"Rohdaten/only21.vcf";
    let mut numberofvariants: i32 = 0;
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
    }
    // ab hier ins neue Modul, kriegt seed, frequenz, genloci und numberofvariants
    let mut sol = HashMap::new();
    let mut beginning: i32 = 0;
    let for_calculation = numberofvariants as f32; 
    println! ("Number of genloci {}", genloci);
    GreedySolver::new().variantselection(genloci, numberofvariants, wanted_freq);
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
    OutputWriter::new().variantselection(sol, genloci, 21);
    let average_time = now.elapsed().as_millis();
    println!("average_time {} ms", average_time);
}






