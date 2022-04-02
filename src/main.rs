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


fn main() {
    // Define problem variables

    //Number of Loci equals the full length of chr21
    let path = &"Rohdaten/only21.vcf";
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

    println!("{:?}", v); //Inhalt des Vectors

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
        },
        Err(msg) => println!("{}", msg),
    }
}


// select the variants and create the output
fn variantselection(solution: Solution){
    println!("in der Funktion");
    println!("Status {:?}", solution.status);
    let mut sorted: Vec<_> = solution.results.iter().collect();
    sorted.sort_by_key(|a| a.0);

    //IndexedReader for fetching
    let third_path = &"Rohdaten/only21.vcf.gz";
    let mut vcf = IndexedReader::from_path(third_path).expect("Error opening file.");
    //let mut to_fetch_from = IndexedReader::from_path(third_path).expect("Error opening file.");
    
   // creating the output vcf
    let head = createheader();
    let mut output = Writer::from_path("Simulationen/test.vcf",&head, true, Format::Vcf).unwrap();
    //Entries that a written into the output
    let mut entries = output.empty_record();
    //Select the entries from the original vcf
    let mut index = 0;
    for line in vcf.records() {
        //println!("{:?}", sorted[index]);
        let coluum = line.expect("Problem with coluum");
        let genotypes =  coluum.genotypes().expect("Error reading genotypes");
        let allel_var = genotypes.get(0);
        //only heterzygotic Variants where the variant is called in CHM1 and the ref is called in CHM13 
        if allel_var[0] != rust_htslib::bcf::record::GenotypeAllele::Unphased(1) || allel_var[1] != rust_htslib::bcf::record::GenotypeAllele::Phased(0){
            continue;
        }
        if *sorted[index].1 == 1.0 {
            let rid = output.header().name2rid(b"21").unwrap(); //still hardcodes b"21"
            entries.set_rid(Some(rid));
            entries.set_pos(coluum.pos());
            entries.set_alleles(&coluum.alleles()).expect("Failed to set alleles");
            let alleles = &[allel_var[0], allel_var[1]];
            entries.push_genotypes(alleles).unwrap(); 
            output.write(&entries).expect("Problem with writing");           
       }
              index += 1;
              //For now because it wont work with the full number of variants
              if index >= sorted.len(){
                  break;
              }
    }
    
    /*Proof of concept for the writer
    for line in vcf.records(){
        let mut coluum = line.expect("Problem with coluum");
        let rid = output.header().name2rid(b"21").unwrap(); //still hardcodes b"21"
        /*let rid = coluum.rid().expect("Kein rid");*/
        entries.set_rid(Some(rid));
        entries.set_pos(coluum.pos());

        let genotypes =  coluum.genotypes().expect("Error reading genotypes");
        let allel_var = genotypes.get(0);
        let alleles = &[allel_var[0], allel_var[1]];
        entries.push_genotypes(alleles).unwrap(); 
        output.write(&entries).unwrap();
    }*/
}


//creates the Header for the output vcf
fn createheader() -> Header {
    let mut header = Header::new();
    let header_contig_line = r#"##contig=<ID=21,length=10>"#;
    header.push_record(header_contig_line.as_bytes());
    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    header.push_record(header_gt_line.as_bytes());
    header.push_sample("syndip".as_bytes());
    header
}
