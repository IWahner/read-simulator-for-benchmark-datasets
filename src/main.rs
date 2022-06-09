use rust_htslib::bcf::*;

use code::linear_solver::LinearSolver;
use code::greedy_solver::GreedySolver;
use code::output_writer::OutputWriter;
use code::errors::Error;
use std::time::{Instant};
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::Read as OtherRead;
use std::collections::HashMap;


#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct SimData {
    variantfreq: f32,
    rng_seed: u64,
    solver: String,
}


fn main() -> Result<(), Error>{
    //time measurment
    let now = Instant::now();

    //reading from the input yaml
    let filename = &"Rohdaten/test.yml";
    let wanted_freq: f32;
    let seed: u64;
    let solver: String;
    match File::open(filename) {
        Ok(mut file) => {
            let mut content = String::new();
            file.read_to_string(&mut content).unwrap();
            let input_yml: SimData = serde_yaml::from_str(&content).unwrap();
            wanted_freq = input_yml.variantfreq;
            seed = input_yml.rng_seed;
            solver = input_yml.solver;
        }
        Err(_) => {
            return Err(Error::FileNotFound);
        }
    }
    if wanted_freq > 1.0 || wanted_freq < 0.0 {
        return Err(Error::InvalidFrequency);
    }
    println!("{}", wanted_freq);
    //Number of Loci equals the full length of chr21
    let path = &"Rohdaten/filteredonlychm21.vcf";
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
    let second_path = &"Rohdaten/filteredonlychm21.vcf";
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
    println! ("Number of genloci {}", genloci);
    let sol: HashMap<String, f32>;
    if solver.eq("linear") {
         sol = LinearSolver::new().variantselection(genloci, numberofvariants, wanted_freq, seed);
    }
    else if solver.eq("greedy") {
         sol = GreedySolver::new().variantselection(genloci, numberofvariants, wanted_freq);
    }
    else {
        return Err(Error::SolverInputProblem);
    }
    OutputWriter::new().write_output(sol, genloci, 21);
    let average_time = now.elapsed().as_millis();
    println!("average_time {} ms", average_time);
    Ok(())
}






