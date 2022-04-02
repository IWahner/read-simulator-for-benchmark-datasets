use rust_htslib::{bam, bam::Read};


pub struct BamWriter {
}

impl<E> BamWriter{

    pub fn new() -> Self {
        BamWriter {
        }
    }

    pub fn write() -> Result<(), E>{
        let mut bam = bam::IndexedReader::from_path(&"Rohdaten/chr21.bam").unwrap();
        Ok(())
    }
}