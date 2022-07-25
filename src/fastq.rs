use std::io::BufRead;
use anyhow::{Result, anyhow};

use super::fastx::FastxRead;
use super::record::Record;

pub struct FastqBytes<B> {
    buf: B
}

impl <B: BufRead> Iterator for FastqBytes<B> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut bytes = Vec::with_capacity(300);
        let mut null = Vec::with_capacity(5);

        let _marker = match self.buf.read_until(b'@', &mut null) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(1) => {},
            Ok(_) => return Some(Err(anyhow!("Misplaced Fastq Marker Sequence '@'")))
        };
        let id = match self.buf.read_until(b'\n', &mut bytes) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(x) => x
        };
        let seq = match self.buf.read_until(b'\n', &mut bytes) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(x) => x
        };
        let plus = match self.buf.read_until(b'\n', &mut bytes) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(x) => x
        };
        let qual = match self.buf.read_until(b'\n', &mut bytes) {
            Err(why) => return Some(Err(anyhow!(why))),
            Ok(0) => return None,
            Ok(x) => x
        };
        let record = Record::new_fastq(bytes, id, seq, plus, qual);
        Some(Ok(record))
    }

}

pub struct FastqReader <R: BufRead> {
    reader: FastqBytes<R>
}
impl <R: BufRead> FastqReader <R> {

    /// # Example
    /// Creates a new [`FastqReader`] explicitly from an object
    /// which implements [`BufRead`].
    ///
    /// ```
    /// let fastq: &'static [u8] = b">sequence.id\nACGTACGT\n+\n$^$%^AA";
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
        Self { reader: FastqBytes { buf: reader } }
    }

    fn next_buffer(&mut self) -> Result<Option<Record>> {
        let buffer = match self.reader.next() {
            Some(line) => Some(line?),
            None => None
        };
        Ok(buffer)
    }
}

impl <R: BufRead> FastxRead for FastqReader<R> {
    fn next_record(&mut self) -> Result<Option<Record>> {
        let buffer = match self.next_buffer()? {
            Some(fastq) => fastq,
            None => return Ok(None)
        };
        Ok(Some(buffer))
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
        let fastq: &'static [u8] = b"@seq.id\nACGT\n+\n7162\n";
        let mut reader = FastqReader::new(fastq);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), b"seq.id");
        assert_eq!(record.as_ref().unwrap().seq(), b"ACGT");
        assert_eq!(reader.into_iter().count(), 0);
    }

    #[test]
    fn unexpected_chars() {
        let fastq: &'static [u8] = b"@seq.id\nABCD\n+\n7162\n";
        let mut reader = FastqReader::new(fastq);
        let record = reader.next().unwrap();
        assert!(!record.valid())
    }

    #[test]
    fn lower_to_upper() {
        let fastq: &'static [u8] = b"@seq.id\nacgt\n+\n7162\n";
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
