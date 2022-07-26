/// An instance of a Fastx Record.
/// This is a two attribute object containing the sequence
/// ID and the Sequence.
#[derive(Debug)]
pub struct Record {
    data: Vec<u8>,
    id: usize,
    seq: usize,
    plus: Option<usize>,
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
            plus: None,
            qual: None
        }
    }

    /// # Usage
    /// Creates a new instance of a `[Record]` from a preinitialized
    /// `[Vec<u8>]` with the `id` and `seq` endpoints calculated. These
    /// endpoints are inclusive of the '\n' terminator and the data is
    /// expected to exclude the prefix '>' marker.
    /// ```
    /// let data = b"seq.0\nACGT\n".to_vec();
    /// let id = 6;
    /// let seq = 5;
    /// let fasta = fxread::Record::new_fasta(data, id, seq);
    /// assert_eq!(fasta.id(), b"seq.0");
    /// assert_eq!(fasta.seq(), b"ACGT");
    /// ```
    pub fn new_fasta(data: Vec<u8>, id: usize, seq: usize) -> Self {
        Self {
            data, id, seq, 
            plus: None, qual: None
        }
    }

    /// # Usage
    /// Creates a new instance of a `[Record]` from a preinitialized
    /// `[Vec<u8>]` with the `id`, `seq`, `plus`, and `qual` endpoints calculated. 
    /// These endpoints are inclusive of the '\n' terminator and the data is
    /// expected to exclude the prefix '@' marker.
    /// ```
    /// let data = b"seq.0\nACGT\n+\n1234\n".to_vec();
    /// let id = 6;
    /// let seq = 5;
    /// let plus = 2;
    /// let qual = 5;
    /// let fasta = fxread::Record::new_fastq(data, id, seq, plus, qual);
    /// assert_eq!(fasta.id(), b"seq.0");
    /// assert_eq!(fasta.seq(), b"ACGT");
    /// assert_eq!(fasta.plus().unwrap(), b"+");
    /// assert_eq!(fasta.qual().unwrap(), b"1234");
    /// ```
    pub fn new_fastq(data: Vec<u8>, id: usize, seq: usize, plus:usize, qual: usize) -> Self {
        Self {
            data, id, seq, 
            plus: Some(plus), qual: Some(qual)
        }
    }

    /// Checks if `[Record]` is empty
    pub fn empty(&self) -> bool {
        (self.id == 0) | (self.seq == 0)
    }

    /// Returns a reference of the sequence ID
    pub fn id(&self) -> &[u8] {
        &self.data[..self.id-1]
    }

    /// Returns a reference of the sequence 
    pub fn seq(&self) -> &[u8] {
        &self.data[self.id..self.id+self.seq-1]
    }

    /// Returns a reference of the '+' region of a fastq
    pub fn plus(&self) -> Option<&[u8]> {
        match self.plus {
            Some(plus) => Some(&self.data[self.id+self.seq..self.id+self.seq+plus-1]),
            None => None
        }
    }

    /// Returns a reference of the sequence 
    pub fn qual(&self) -> Option<&[u8]> {
        let plus = match self.plus {
            Some(plus) => plus,
            None => return None
        };
        match self.qual {
            Some(qual) => Some(&self.data[self.id+self.seq+plus..self.id+self.seq+plus+qual-1]),
            None => None
        }
    }

    /// Returns a reference to the raw data underlying the record
    pub fn data(&self) -> &[u8] {
        &self.data
    }

    /// Validates that the record is not partially constructed
    /// or composed of unexpected nucleotides
    pub fn valid(&self) -> bool {
        if self.empty() { false }
        else { self.valid_sequence() }
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
        self.seq().iter().all(|b| matches!(b, b'A'|b'a'|b'C'|b'c'|b'G'|b'g'|b'T'|b't'|b'N'|b'n'|b'U'|b'u'))
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod test {
    use super::Record;

    fn gen_valid_fasta() -> (Vec<u8>, usize, usize) {
        (b"seq.0\nACGT\n".to_vec(), 6, 5)
    }

    fn gen_invalid_fasta() -> (Vec<u8>, usize, usize) {
        (b"seq.0\nABCD\n".to_vec(), 6, 5)
    }

    fn gen_valid_fasta_lower() -> (Vec<u8>, usize, usize) {
        (b"seq.0\nacgt\n".to_vec(), 6, 5)
    }

    fn gen_invalid_fasta_lower() -> (Vec<u8>, usize, usize) {
        (b"seq.0\nabcd\n".to_vec(), 6, 5)
    }

    fn gen_valid_fasta_rev() -> (Vec<u8>, usize, usize) {
        (b"seq.0\nATCGGCTA\n".to_vec(), 6, 9)
    }

    fn gen_valid_fasta_rev_lower() -> (Vec<u8>, usize, usize) {
        (b"seq.0\natcggcta\n".to_vec(), 6, 9)
    }

    fn gen_valid_fastq() -> (Vec<u8>, usize, usize, usize, usize) {
        (b"seq.0\nACGT\n+\n1234\n".to_vec(), 6, 5, 2, 5)
    }

    fn gen_invalid_fastq() -> (Vec<u8>, usize, usize, usize, usize) {
        (b"seq.0\nABCD\n+\n1234\n".to_vec(), 6, 5, 2, 5)
    }

    #[test]
    fn create() {
        let record = Record::new();
        assert!(record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn valid_fasta() {
        let (fasta, id, seq) = gen_valid_fasta();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.id(), b"seq.0");
        assert_eq!(record.seq(), b"ACGT");
    }

    #[test]
    fn invalid_fasta() {
        let (fasta, id, seq) = gen_invalid_fasta();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn valid_fastq() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.id(), b"seq.0");
        assert_eq!(record.seq(), b"ACGT");
        assert_eq!(record.plus().unwrap(), b"+");
        assert_eq!(record.qual().unwrap(), b"1234");
    }

    #[test]
    fn invalid_fastq() {
        let (fasta, id, seq, plus, qual) = gen_invalid_fastq();
        let record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(!record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn valid_lowercase() {
        let (fasta, id, seq) = gen_valid_fasta_lower();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.id(), b"seq.0");
        assert_eq!(record.seq(), b"acgt");
    }

    #[test]
    fn invalid_lowercase() {
        let (fasta, id, seq) = gen_invalid_fasta_lower();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn upper_conversion() {
        let (fasta, id, seq) = gen_valid_fasta_lower();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq_upper(), b"ACGT");
    }

    #[test]
    fn reverse_complement() {
        let (fasta, id, seq) = gen_valid_fasta_rev();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq_rev_comp(), b"TAGCCGAT");
    }

    #[test]
    fn reverse_complement_lower() {
        let (fasta, id, seq) = gen_valid_fasta_rev_lower();
        let record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq_rev_comp(), b"tagccgat");
    }
}
