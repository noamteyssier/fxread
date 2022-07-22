use std::io::BufRead;
use anyhow::{anyhow, Result};
use bstr::{ByteSlice, io::{ByteLines, BufReadExt}};

use super::fastx::FastxRead;
use super::record::Record;

/// A Fastq Reader implementation.
pub struct FastqReader <R: BufRead> {
    reader: ByteLines<R>,
}
impl <R: BufRead> FastqReader <R> {

    /// # Example
    /// Creates a new [`FastqReader`] explicitly from an object
    /// which implements [`BufRead`].
    ///
    /// ```
    /// let fastq: &'static [u8] = b">sequence.id\nACGTACGT\n+\n$^$%^2@@";
    /// let reader = fxread::FastqReader::new(fastq);
    /// ```
    ///
    /// or from a file reader.
    /// ```
    /// let file = std::fs::File::open("example/sequences.fq").unwrap();
    /// let buffer = std::io::BufReader::new(file);
    /// let reader = fxread::FastqReader::new(buffer);
    /// ```
    pub fn new(reader: R) -> Self {
        Self { reader: reader.byte_lines() }
    }

    fn parse_header(token: &[u8]) -> Result<&[u8]> {
        match token[0] {
            b'@' => Ok(&token[1..]),
            _ => Err(anyhow!("Header does not begin with '@' "))
        }
    }
}

impl <R: BufRead> FastxRead for FastqReader<R> {
    fn next_record(&mut self) -> Result<Option<Record>> {
        let mut record = Record::new();

        for idx in 0..4 {
            let buffer = match self.reader.next() {
                Some(line) => line,
                None => return Ok(None)
            }?;
            match idx {
                0 => {
                    record.set_id(Self::parse_header(&buffer)?)
                },
                1 => {
                    record.set_seq(&buffer)
                },
                _ => {}
            }
        }

        if record.empty() {
            Ok(None)
        }
        else {
            Ok(Some(record))
        }
    }
}

impl <R: BufRead> Iterator for FastqReader <R> {

    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(r) => r,
            Err(why) => panic!("{}", why)
        }
    }

}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use flate2::read::MultiGzDecoder;
    use super::FastqReader;
    
    #[test]
    fn read_string() {
        let fastq: &'static [u8] = b"@seq.id\nACGT\n+\n7162";
        let mut reader = FastqReader::new(fastq);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.id");
        assert_eq!(record.as_ref().unwrap().seq(), b"ACGT");
        assert_eq!(reader.into_iter().count(), 0);
    }

    #[test]
    fn unexpected_chars() {
        let fastq: &'static [u8] = b"@seq.id\nABCD\n+\n7162";
        let mut reader = FastqReader::new(fastq);
        let record = reader.next().unwrap();
        assert!(!record.valid())
    }

    #[test]
    fn lower_to_upper() {
        let fastq: &'static [u8] = b"@seq.id\nacgt\n+\n7162";
        let mut reader = FastqReader::new(fastq);
        let record = reader.next().unwrap();
        assert_eq!(record.seq_upper(), b"ACGT");
    }

    #[test]
    fn read_plaintext() {
        let file = File::open("example/sequences.fq").unwrap();
        let buffer = BufReader::new(file);
        let mut reader = FastqReader::new(buffer);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.0");
        assert_eq!(record.as_ref().unwrap().seq(), b"TAGTGCTTTCGATGGAACTGGACCGAGAATTCTATCGCAAATGGAACCGGAGTGACGGTGTTTCTAGACGCTCCTCACAA");
        assert_eq!(reader.into_iter().count(), 9);
    }

    #[test]
    fn read_gzip() {
        let file = File::open("example/sequences.fq.gz").unwrap();
        let gzip = MultiGzDecoder::new(file);
        let buffer = BufReader::new(gzip);
        let mut reader = FastqReader::new(buffer);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.0");
        assert_eq!(record.as_ref().unwrap().seq(), b"TAGTGCTTTCGATGGAACTGGACCGAGAATTCTATCGCAAATGGAACCGGAGTGACGGTGTTTCTAGACGCTCCTCACAA");
        assert_eq!(reader.into_iter().count(), 9);
    }
}
