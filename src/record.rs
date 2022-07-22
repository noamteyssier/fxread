use anyhow::{anyhow, Result};
use memchr::memchr;

/// An instance of a Fastx Record.
/// This is a two attribute object containing the sequence
/// ID and the Sequence.
#[derive(Debug)]
pub struct Record {
    data: Vec<u8>,
    id: usize,
    seq: usize,
    qual: Option<usize>
}
impl Record {

    /// # Usage
    /// Creates a new instance of a `[Record]`
    /// ```
    /// let record = fxread::Record::new();
    /// assert!(record.empty());
    /// ```
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            id: 0,
            seq: 0,
            qual: None
        }
    }

    /// Initializes a record from fasta data
    pub fn from_fasta(data: Vec<u8>) -> Result<Self> {
        let id = match memchr(b'\n', &data) {
            Some(index) => index,
            None => return Err(anyhow!("Unexpected Termination in header"))
        };
        let seq = match memchr(b'\n', &data[id+1..]) {
            Some(index) => index,
            None => return Err(anyhow!("Unexpected Termination in sequence"))
        };
        Ok(
            Self {
                data,
                id,
                seq,
                qual: None
            }
        )
    }

    /// Initializes a record from fastq data
    pub fn from_fastq(data: Vec<u8>) -> Result<Self> {

        let id = match memchr(b'\n', &data) {
            Some(index) => index,
            None => return Err(anyhow!("Unexpected Termination in Header"))
        };

        let seq = match memchr(b'\n', &data[id+1..]) {
            Some(index) => index,
            None => return Err(anyhow!("Unexpected Termination in Sequence"))
        };

        let plus = match memchr(b'\n', &data[id+seq+2..]) {
            Some(index) => index,
            None => return Err(anyhow!("Unexpected Termination in Plus"))
        };

        let qual = match memchr(b'\n', &data[id+seq+plus+3..]) {
            Some(index) => index,
            None => return Err(anyhow!("Unexpected Termination in Qual"))
        };

        Ok(
            Self {
                data,
                id,
                seq,
                qual: Some(qual)
            }
        )
    }

    /// Checks if `[Record]` is empty
    pub fn empty(&self) -> bool {
        (self.id == 0) | (self.seq == 0)
    }

    /// Returns a reference of the sequence ID
    pub fn id(&self) -> &[u8] {
        &self.data[..self.id]
    }

    /// Returns a reference of the sequence 
    pub fn seq(&self) -> &[u8] {
        &self.data[self.id+1..self.id+1+self.seq]
    }

    /// Returns a reference of the sequence 
    pub fn qual(&self) -> Option<&[u8]> {
        match self.qual {
            Some(qual) => Some(&self.data[self.id+self.seq+4..self.id+self.seq+4+qual]),
            None => None
        }
    }

    /// Validates that the record is not partially constructed
    /// or composed of unexpected nucleotides
    pub fn valid(&self) -> bool {
        !self.empty() & self.valid_sequence()
    }

    /// Converts the sequence to uppercase
    pub fn seq_upper(&self) -> Vec<u8> {
        self.seq()
            .iter()
            .map(|c| if c & b' ' == 0 { *c } else { c ^ b' ' })
            .collect()
    }

    /// Reverse Complements the sequence
    pub fn seq_rev_comp(&self) -> Vec<u8> {
        self.seq()
            .iter()
            .rev()
            .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
            .collect()
    }

    /// Validates whether sequence is composed
    /// of valid nucleotides
    fn valid_sequence(&self) -> bool {
        self.seq().iter().all(|b| match b {
            b'A'|b'a'|b'C'|b'c'|b'G'|b'g'|b'T'|b't'|b'N'|b'n'|b'U'|b'u' => true,
            _ => false
        })
    }

}

#[cfg(test)]
mod test {
    //use super::Record;

    //#[test]
    //fn create() {
        //let record = Record::new();
        //assert!(record.empty());
        //assert!(!record.valid());
    //}

    //#[test]
    //fn create_partial_id() {
        //let mut record = Record::new();
        //record.set_id(b"some_id".to_vec());
        //assert!(record.empty());
        //assert!(!record.valid());
    //}

    //#[test]
    //fn create_partial_seq() {
        //let mut record = Record::new();
        //record.set_seq(b"ACGT".to_vec());
        //assert!(record.empty());
        //assert!(!record.valid());
    //}

    //#[test]
    //fn valid() {
        //let mut record = Record::new();
        //record.set_id(b"some_id".to_vec());
        //record.set_seq(b"ACGT".to_vec());
        //assert!(!record.empty());
        //assert!(record.valid());
    //}

    //#[test]
    //fn invalid() {
        //let mut record = Record::new();
        //record.set_id(b"some_id".to_vec());
        //record.set_seq(b"BCGT".to_vec());
        //assert!(!record.empty());
        //assert!(!record.valid());
    //}

    //#[test]
    //fn valid_lowercase() {
        //let mut record = Record::new();
        //record.set_id(b"some_id".to_vec());
        //record.set_seq(b"acgt".to_vec());
        //assert!(!record.empty());
        //assert!(record.valid());
    //}

    //#[test]
    //fn invalid_lowercase() {
        //let mut record = Record::new();
        //record.set_id(b"some_id".to_vec());
        //record.set_seq(b"bcgt".to_vec());
        //assert!(!record.empty());
        //assert!(!record.valid());
    //}

    //#[test]
    //fn upper_conversion() {
        //let mut record = Record::new();
        //record.set_seq(b"acgt".to_vec());
        //assert_eq!(record.seq_upper(), b"ACGT");
    //}

    //#[test]
    //fn reverse_complement() {
        //let mut record = Record::new();
        //record.set_seq(b"ACGTA".to_vec());
        //assert_eq!(record.seq_rev_comp(), b"TACGT");
    //}

    //#[test]
    //fn lower_reverse_complement() {
        //let mut record = Record::new();
        //record.set_seq(b"acgta".to_vec());
        //assert_eq!(record.seq_rev_comp(), b"tacgt");
    //}
}
