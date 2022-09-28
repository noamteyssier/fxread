use super::fastx::FastxRead;
use super::record::Record;
use anyhow::{anyhow, Result};
use std::io::BufRead;

/// Struct to handle the Byte Reading for Fasta Formatted Files.
/// Heavily inspired from bstr `ByteRecord`.
pub struct FastaBytes<B> {
    buf: B,
}

impl<B: BufRead> Iterator for FastaBytes<B> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut bytes = Vec::with_capacity(300);
        let mut null = Vec::with_capacity(5);

        match self.buf.read_until(b'>', &mut null) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(1) => {}
            Ok(_) => return Some(Err(anyhow!("Misplaced Fasta Marker Sequence '>'"))),
        };
        let id = match self.buf.read_until(b'\n', &mut bytes) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(x) => x,
        };
        let seq = match self.buf.read_until(b'\n', &mut bytes) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(x) => x,
        };
        let record = Record::new_fasta(bytes, id, seq);
        Some(Ok(record))
    }
}

/// A Fasta Reader implementation.
pub struct FastaReader<R: BufRead> {
    reader: FastaBytes<R>,
}
impl<R: BufRead> FastaReader<R> {
    /// # Example
    /// Creates a new [`FastaReader`] explicitly from an object
    /// which implements [`BufRead`].
    ///
    /// ```
    /// let fasta: &'static [u8] = b">sequence.id\nACGTACGT\n";
    /// let reader = fxread::FastaReader::new(fasta);
    /// ```
    ///
    /// or from a file reader.
    /// ```
    /// let file = std::fs::File::open("example/sequences.fa").unwrap();
    /// let buffer = std::io::BufReader::new(file);
    /// let reader = fxread::FastaReader::new(buffer);
    /// ```
    pub fn new(reader: R) -> Self {
        Self {
            reader: FastaBytes { buf: reader },
        }
    }

    fn next_buffer(&mut self) -> Result<Option<Record>> {
        let buffer = match self.reader.next() {
            Some(line) => Some(line?),
            None => None,
        };
        Ok(buffer)
    }
}

impl<R: BufRead> FastxRead for FastaReader<R> {
    fn next_record(&mut self) -> Result<Option<Record>> {
        let buffer = match self.next_buffer()? {
            Some(fastq) => fastq,
            None => return Ok(None),
        };
        Ok(Some(buffer))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(r) => r,
            Err(why) => panic!("{}", why),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::FastaReader;
    use flate2::read::MultiGzDecoder;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn read_string() {
        let fasta: &'static [u8] = b">seq.id\nACGT\n";
        let mut reader = FastaReader::new(fasta);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.id");
        assert_eq!(record.as_ref().unwrap().seq(), b"ACGT");
        assert_eq!(reader.into_iter().count(), 0);
    }

    #[test]
    fn unexpected_chars() {
        let fasta: &'static [u8] = b">seq.id\nABCD\n";
        let mut reader = FastaReader::new(fasta);
        let record = reader.next().unwrap();
        assert!(!record.valid())
    }

    #[test]
    fn lower_to_upper() {
        let fasta: &'static [u8] = b">seq.id\nacgt\n";
        let mut reader = FastaReader::new(fasta);
        let record = reader.next().unwrap();
        assert_eq!(record.seq_upper(), b"ACGT");
    }

    #[test]
    fn read_plaintext() {
        let file = File::open("example/sequences.fa").unwrap();
        let buffer = BufReader::new(file);
        let mut reader = FastaReader::new(buffer);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.0");
        assert_eq!(
            record.as_ref().unwrap().seq(),
            b"TAGTGCTTTCGATGGAACTGGACCGAGAATTCTATCGCAAATGGAACCGGAGTGACGGTGTTTCTAGACGCTCCTCACAA"
        );
        assert_eq!(reader.into_iter().count(), 9);
    }

    #[test]
    fn read_gzip() {
        let file = File::open("example/sequences.fa.gz").unwrap();
        let gzip = MultiGzDecoder::new(file);
        let buffer = BufReader::new(gzip);
        let mut reader = FastaReader::new(buffer);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.0");
        assert_eq!(
            record.as_ref().unwrap().seq(),
            b"TAGTGCTTTCGATGGAACTGGACCGAGAATTCTATCGCAAATGGAACCGGAGTGACGGTGTTTCTAGACGCTCCTCACAA"
        );
        assert_eq!(reader.into_iter().count(), 9);
    }
}
