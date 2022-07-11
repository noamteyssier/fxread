use std::{io::{BufReader, BufRead}, fs::File};
use flate2::read::MultiGzDecoder;
use anyhow::Result;

mod record;
mod fastx;
mod fastq;
mod fasta;

use fastx::FastxRead;
use fasta::FastaReader;
use fastq::FastqReader;
use record::Record;


fn initialize_generic_buffer(path: &str, is_gzip: bool) -> Result<Box<dyn BufRead>> {
    match is_gzip {
        true => {
            let file = File::open(path)?;
            let gzip = MultiGzDecoder::new(file);
            let buffer = BufReader::new(gzip);
            Ok(Box::new(buffer))
        },
        false => {
            let file = File::open(path)?;
            let buffer = BufReader::new(file);
            Ok(Box::new(buffer))
        }
    }
}

fn initialize_generic_reader(buffer: Box<dyn BufRead>, is_fasta: bool) -> Box<dyn FastxRead<Item = Record>> {
    match is_fasta {
        true => Box::new(FastaReader::new(buffer)),
        false => Box::new(FastqReader::new(buffer))
    }
}


fn initialize_reader(path: &str) -> Result<Box<dyn FastxRead<Item = Record>>>{
    let is_gzip = path.ends_with(".gz");
    let is_fasta = path.contains(".fa") | path.contains(".fasta");
    let buffer = initialize_generic_buffer(path, is_gzip)?;
    let reader = initialize_generic_reader(buffer, is_fasta);
    Ok(reader)
}


fn run_reader(reader: Box<dyn FastxRead<Item = Record>>) {
    reader
        .into_iter()
        .for_each(|x| println!("{:?}", x));
}


fn main() -> Result<()> {
    let path = "sequences.fq.gz";
    let reader = initialize_reader(path)?;
    run_reader(reader);

    Ok(())
}
