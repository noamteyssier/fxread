use anyhow::Result;
use std::{
    convert::AsRef,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use super::{FastaReader, FastqReader, FastxRead, Record};

const BUFFER_SIZE: usize = 4096 * 68;

fn initialize_generic_buffer<P>(path: P) -> Result<Box<BufReader<Box<dyn std::io::Read>>>>
where
    P: AsRef<Path>,
{
    Ok(Box::new(std::io::BufReader::new(
        niffler::get_reader(Box::new(File::open(path)?))?.0,
    )))
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
pub fn initialize_reader<P>(path: P) -> Result<Box<dyn FastxRead<Item = Record>>>
where
    P: AsRef<Path>,
{
    let mut buffer = initialize_generic_buffer(path)?;
    buffer.fill_buf()?;
    if buffer.buffer().is_empty() {
        return Err(anyhow::anyhow!("No data in input file"));
    }
    match buffer.buffer()[0] {
        b'>' => Ok(initialize_generic_reader(buffer, true)),
        b'@' => Ok(initialize_generic_reader(buffer, false)),
        _ => Err(anyhow::anyhow!("Unrecognized file format")),
    }
}

/// Initializes a reader from stdin. This is useful for piping
/// in data from other programs.
///
/// ## From Stdin
/// This example shows the creation of a reader from stdin.
/// ```
/// use std::io::stdin;
/// use fxread::initialize_stdin_reader;
///
/// let input = stdin().lock();
/// let reader = initialize_stdin_reader(input);
/// ```
pub fn initialize_stdin_reader<R: BufRead + 'static>(
    reader: R,
) -> Result<Box<dyn FastxRead<Item = Record>>> {
    let mut buffer = BufReader::with_capacity(BUFFER_SIZE, reader);
    buffer.fill_buf()?;
    if buffer.buffer().is_empty() {
        return Err(anyhow::anyhow!("No data in stdin"));
    }
    match buffer.buffer()[0] {
        b'>' => Ok(initialize_generic_reader(Box::new(buffer), true)),
        b'@' => Ok(initialize_generic_reader(Box::new(buffer), false)),
        _ => Err(anyhow::anyhow!("Unrecognized file format")),
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use std::io::Cursor;

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
    fn assign_bz2fasta() {
        let path = "example/sequences.fa.bz2";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_xzfasta() {
        let path = "example/sequences.fa.xz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_zstfasta() {
        let path = "example/sequences.fa.zst";
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

    #[test]
    fn assign_bz2fastq() {
        let path = "example/sequences.fq.bz2";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_xzfastq() {
        let path = "example/sequences.fq.xz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_zstfastq() {
        let path = "example/sequences.fq.zst";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader.into_iter().map(|x| assert!(!x.empty())).count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn assign_fa_stdin() {
        let example_fa = ">test\nACGT\n>test2\nACGT\n";
        let cursor = Cursor::new(example_fa);
        let reader = initialize_stdin_reader(cursor).expect("invalid path");
        let num_records = reader.count();
        assert_eq!(num_records, 2);
    }

    #[test]
    fn assign_fq_stdin() {
        let example_fq = "@test\nACGT\n+\n!!!!\n@test2\nACGT\n+\n!!!!\n";
        let cursor = Cursor::new(example_fq);
        let reader = initialize_stdin_reader(cursor).expect("invalid path");
        let num_records = reader.count();
        assert_eq!(num_records, 2);
    }

    #[test]
    fn assign_malformed_stdin() {
        let example_malformed = "test\nACGT\n+\n!!!!\n@test2\nACGT\n+\n!!!!\n";
        let cursor = Cursor::new(example_malformed);
        let reader = initialize_stdin_reader(cursor);
        assert!(reader.is_err());
    }
}
