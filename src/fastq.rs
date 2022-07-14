use std::io::BufRead;
use anyhow::Result;
use super::fastx::FastxRead;
use super::record::Record;

/// A Fastq Reader implementation.
pub struct FastqReader <R: BufRead> {
    reader: R,
    buffer: String
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
        Self { 
            reader,
            buffer: String::new()
        }
    }

    fn next_line(&mut self) -> Result<bool> {
        Ok(self.reader.read_line(&mut self.buffer)? > 0)
    }

    fn strip_header(&self, token: &str) -> String {
        token
            .trim_start_matches('@')
            .trim_end()
            .to_string()
    }

    fn strip_sequence(&self, token: &str) -> String {
        token
            .trim_end()
            .to_string()
    }
}

impl <R: BufRead> FastxRead for FastqReader<R> {
    fn next_record(&mut self) -> Result<Option<Record>> {
        let mut record = Record::new();

        for idx in 0..4 {
            self.buffer.clear();
            if !self.next_line()? { break }
            match idx {
                0 => {
                    record.set_id(self.strip_header(&self.buffer))
                },
                1 => {
                    record.set_seq(self.strip_sequence(&self.buffer))
                },
                _ => { }
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
            Err(_) => panic!("Unexpected file end")
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
        let fasta: &'static [u8] = b"@seq.id\nACGT\n+\n7162";
        let mut reader = FastqReader::new(fasta);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), "seq.id");
        assert_eq!(record.as_ref().unwrap().seq(), "ACGT");
        assert_eq!(reader.into_iter().count(), 0);
    }

    #[test]
    fn read_plaintext() {
        let file = File::open("example/sequences.fq").unwrap();
        let buffer = BufReader::new(file);
        let mut reader = FastqReader::new(buffer);
        let record = reader.next();
        assert!(record.as_ref().is_some());
        assert_eq!(record.as_ref().unwrap().id(), "seq.0");
        assert_eq!(record.as_ref().unwrap().seq(), "TAGTGCTTTCGATGGAACTGGACCGAGAATTCTATCGCAAATGGAACCGGAGTGACGGTGTTTCTAGACGCTCCTCACAA");
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
        assert_eq!(record.as_ref().unwrap().id(), "seq.0");
        assert_eq!(record.as_ref().unwrap().seq(), "TAGTGCTTTCGATGGAACTGGACCGAGAATTCTATCGCAAATGGAACCGGAGTGACGGTGTTTCTAGACGCTCCTCACAA");
        assert_eq!(reader.into_iter().count(), 9);
    }
}
