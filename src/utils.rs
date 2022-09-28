use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use super::{FastaReader, FastqReader, FastxRead, Record};

const BUFFER_SIZE: usize = 4096 * 68;

/// Handles the algebraic enumerations of possible expected file types.
/// This is either a Fasta/Fastq and whether or not it is Gzipped.
#[derive(Debug)]
enum FileType {
    Fasta,
    Fastq,
    FastaGz,
    FastqGz,
}
impl FileType {
    /// Handles creation of enum from path name
    pub fn from_path(path: &str) -> Option<Self> {
        if path.ends_with(".fastq.gz") | path.ends_with(".fq.gz") {
            Some(Self::FastqGz)
        } else if path.ends_with(".fastq") | path.ends_with(".fq") {
            Some(Self::Fastq)
        } else if path.ends_with(".fasta.gz") | path.ends_with(".fa.gz") {
            Some(Self::FastaGz)
        } else if path.ends_with(".fasta") | path.ends_with(".fa") {
            Some(Self::Fasta)
        } else {
            None
        }
    }
}

fn initialize_generic_buffer(path: &str, is_gzip: bool) -> Result<Box<dyn BufRead>> {
    if is_gzip {
        let file = File::open(path)?;
        let gzip = MultiGzDecoder::new(file);
        let buffer = BufReader::with_capacity(BUFFER_SIZE, gzip);
        Ok(Box::new(buffer))
    } else {
        let file = File::open(path)?;
        let buffer = BufReader::with_capacity(BUFFER_SIZE, file);
        Ok(Box::new(buffer))
    }
}

fn initialize_generic_reader(
    buffer: Box<dyn BufRead>,
    is_fasta: bool,
) -> Box<dyn FastxRead<Item = Record>> {
    if is_fasta {
        Box::new(FastaReader::new(buffer))
    } else {
        Box::new(FastqReader::new(buffer))
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
///     .for_each(|record| println!("{:?}", record));
/// ```
pub fn initialize_reader(path: &str) -> Result<Box<dyn FastxRead<Item = Record>>> {
    let (is_gzip, is_fasta) = match FileType::from_path(path) {
        Some(s) => match s {
            FileType::Fasta => Ok((false, true)),
            FileType::Fastq => Ok((false, false)),
            FileType::FastaGz => Ok((true, true)),
            FileType::FastqGz => Ok((true, false)),
        },
        None => Err(anyhow::anyhow!(
            "Cannot accurately parse filetype from path: {}",
            path
        )),
    }?;
    let buffer = initialize_generic_buffer(path, is_gzip)?;
    let reader = initialize_generic_reader(buffer, is_fasta);
    Ok(reader)
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn assignment_fastq_short() {
        let path = "sequence.fq";
        match FileType::from_path(path) {
            Some(FileType::Fastq) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_fastq_long() {
        let path = "sequence.fastq";
        match FileType::from_path(path) {
            Some(FileType::Fastq) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_fasta_short() {
        let path = "sequence.fa";
        match FileType::from_path(path) {
            Some(FileType::Fasta) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_fasta_long() {
        let path = "sequence.fasta";
        match FileType::from_path(path) {
            Some(FileType::Fasta) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_gz_fastq_short() {
        let path = "sequence.fq.gz";
        match FileType::from_path(path) {
            Some(FileType::FastqGz) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_gz_fastq_long() {
        let path = "sequence.fastq.gz";
        match FileType::from_path(path) {
            Some(FileType::FastqGz) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_gz_fasta_short() {
        let path = "sequence.fa.gz";
        match FileType::from_path(path) {
            Some(FileType::FastaGz) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assignment_gz_fasta_long() {
        let path = "sequence.fasta.gz";
        match FileType::from_path(path) {
            Some(FileType::FastaGz) => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn assign_fasta() {
        let path = "example/sequences.fa";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_gzfasta() {
        let path = "example/sequences.fa.gz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_fastq() {
        let path = "example/sequences.fq";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_gzfastq() {
        let path = "example/sequences.fq.gz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }
}
