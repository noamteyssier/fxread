use std::{io::{BufRead, BufReader}, fs::File};
use flate2::read::MultiGzDecoder;
use anyhow::Result;

use super::{
    FastxRead,
    FastaReader,
    FastqReader,
    Record
};

fn initialize_generic_buffer(
        path: &str, 
        is_gzip: bool) -> Result<Box<dyn BufRead>> 
{
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

fn initialize_generic_reader(
        buffer: Box<dyn BufRead>, 
        is_fasta: bool) -> Box<dyn FastxRead<Item = Record>> 
{
    match is_fasta {
        true => Box::new(FastaReader::new(buffer)),
        false => Box::new(FastqReader::new(buffer))
    }
}

pub fn initialize_reader(path: &str) -> Result<Box<dyn FastxRead<Item = Record>>>{
    let is_gzip = path.ends_with(".gz");
    let is_fasta = path.contains(".fa") | path.contains(".fasta");
    let buffer = initialize_generic_buffer(path, is_gzip)?;
    let reader = initialize_generic_reader(buffer, is_fasta);
    Ok(reader)
}


#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn assign_fasta() {
        let path = "example/sequences.fa";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_gzfasta() {
        let path = "example/sequences.fa.gz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_fastq() {
        let path = "example/sequences.fq";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_gzfastq() {
        let path = "example/sequences.fq.gz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }
}
