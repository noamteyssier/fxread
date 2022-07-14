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
        }
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

/// # Initializing a reader dependent on the file path extensions.
/// ## Recognized Extensions
/// This recognizes `FASTA` formats from `*.fa` and `*.fasta` and
/// `FASTQ` formats from `*.fq` and `*.fastq`. These will then be
/// opened with gzip buffers or standard depending on if the 
/// `*.gz` extension is found at the end of the pathname.
/// 
/// ## From Fasta
/// This example shows the creation of a reader from a fasta formatted
/// plaintext file.
/// ```
/// use fxread::initialize_reader;
/// let path = "example/sequences.fa";
/// let reader = initialize_reader(path).unwrap();
/// reader
///     .into_iter()
///     .for_each(|record| println!("{:?}", record));
/// ```
///
/// ## From Gzip Fasta
/// A common problem for the above however is that your file
/// may actually be gzipped. Here you can use the same function
/// and it will handle the initialization of the reader depending
/// on the extension.
/// ```
/// use fxread::initialize_reader;
/// let path = "example/sequences.fa.gz";
/// let reader = initialize_reader(path).unwrap();
/// reader
///     .into_iter()
///     .for_each(|record| println!("{:?}", record));
/// ```
///
/// ## From Fastq
/// Here is another example for a fastq-formatted file that 
/// is plaintext
/// ```
/// use fxread::initialize_reader;
/// let path = "example/sequences.fq";
/// let reader = initialize_reader(path).unwrap();
/// reader
///     .into_iter()
///     .for_each(|record| println!("{:?}", record));
///
/// ```
/// ## From Gzip Fastq
/// Here is another example for a fastq-formatted file that 
/// is gzipped
/// ```
/// use fxread::initialize_reader;
/// let path = "example/sequences.fq.gz";
/// let reader = initialize_reader(path).unwrap();
/// reader
///     .into_iter()
///     .for_each(|record| println!("{:?}", record));
/// ```
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
