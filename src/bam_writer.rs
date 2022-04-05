use rust_htslib::{bam, bam::Read};
use core::fmt::Error;

pub struct BamWriter {
}

impl BamWriter
{

    pub fn new() -> Self {
        BamWriter {
        }
    }


    pub fn write(&self, chr: i32, position: i32) -> Result<(), Error>{
        let mut bam = bam::IndexedReader::from_path(&"Rohdaten/inchr21.bam").unwrap();
        let header = bam::Header::from_template(bam.header());
        let mut out = bam::Writer::from_path(&"Simulationen/test.bam", &header, bam::Format::Bam).unwrap();

        //println!("{:?}", bam.header().target_count());
        //println!("{}", bam.records().count());
        bam.fetch((chr, position-500, position+500)).expect("Fehler!!!"); //for some reason you get reads from chr21, when you use 20 to fetch them
        for read in bam.rc_records() {
            let record = read.unwrap();
            out.write(&record).unwrap();
            //println!("read name: {:?}", read.unwrap().qname());
        }
        Ok(())
    }
}
