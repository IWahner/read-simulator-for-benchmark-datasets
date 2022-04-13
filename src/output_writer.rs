use std::collections::HashMap;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Header;
use rust_htslib::bam::Writer;
use rust_htslib::bam::Read;
use rust_htslib::bam::Format;
use rust_htslib::bcf::IndexedReader as OtherIndexedReader;
use rust_htslib::bcf::Writer as OtherWriter;
use rust_htslib::bcf::Header as OtherHeader;
use rust_htslib::bcf::Read as OtherRead;
use rust_htslib::bcf::Format as OtherFormat;


pub struct OutputWriter {
}

impl Default for OutputWriter {
    fn default() -> Self {
         Self::new()
      }
  }
 
impl OutputWriter{
 
     pub fn new() -> Self {
         OutputWriter {
         }
    }

    // select the variants and creates the vcf and bam output file
    //TODO length and chr as Vec of tuple
    pub fn write_output(&self, solution: HashMap<String, f32>, length: i32, chr: i32){
        //println!("Status {:?}", solution.status);
        let mut sorted: Vec<_> = solution.iter().collect();
        sorted.sort_by_key(|a| a.0);

        //IndexedReader for fetching
        let third_path = &"Rohdaten/filtered.vcf.gz";
        let mut vcf = OtherIndexedReader::from_path(third_path).expect("Error opening file.");
        
    // creating the output vcf and bam
        let head_vcf = self.createheader(length, chr);
        let mut output = OtherWriter::from_path("Simulationen/test.vcf",&head_vcf, true, OtherFormat::Vcf).unwrap();
    
        let mut bam_read = IndexedReader::from_path(&"Rohdaten/inchr21.bam").unwrap();
        let head_bam = Header::from_template(bam_read.header());
        let mut out = Writer::from_path(&"Simulationen/test.bam", &head_bam, Format::Bam).unwrap();

        
        //Entries that a written into the output
        let mut entries = output.empty_record();
        //Select the entries from the original vcf
        let mut index = 0;
        //for some reason you get reads from chr21, when you use 20 to fetch them
        let fetch_chr = chr-1;
        println!("{:?}", sorted[index]);
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
                //println!("{:?}", sorted[index]);
                let rid = output.header().name2rid(b"21").unwrap(); //still hardcodes b"21"
                entries.set_rid(Some(rid));
                let position = coluum.pos();
                entries.set_pos(position);
                entries.set_alleles(&coluum.alleles()).expect("Failed to set alleles");
                let alleles = &[allel_var[0], allel_var[1]];
                entries.push_genotypes(alleles).unwrap(); 
                output.write(&entries).expect("Problem with writing");
                //write the bam file
                bam_read.fetch((fetch_chr, position-500, position+500)).expect("Fehler!!!"); 
                for read in bam_read.rc_records() {
                    let record = read.unwrap();
                    out.write(&record).unwrap();         
                }
            }
            index += 1;
                //For now because it wont work with the full number of variants
            if index >= solution.len(){
                break;
            }
        }
    }
    
       
    //creates the Header for the output vcf
    pub fn createheader(&self, length: i32, chr: i32) -> OtherHeader {
        let mut header = OtherHeader::new();
        let header_contig_line = format!("##contig=<ID={},length={}>", chr, length); 
        println! ("{}", header_contig_line);
        header.push_record(header_contig_line.as_bytes());
        let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
        header.push_record(header_gt_line.as_bytes());
        header.push_sample("syndip".as_bytes());
        header
    }

}    
