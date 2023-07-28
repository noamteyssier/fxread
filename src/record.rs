use std::ops::Range;

/// An instance of a Fastx Record.
/// This is a two attribute object containing the sequence
/// ID and the Sequence.
#[derive(Debug)]
pub struct Record {
    data: Vec<u8>,
    id: usize,
    seq: usize,
    plus: Option<usize>,
    qual: Option<usize>,
}
impl Record {
    /// # Usage
    /// Creates a new instance of a `[Record]`
    /// ```
    /// let record = fxread::Record::new();
    /// assert!(record.empty());
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            id: 0,
            seq: 0,
            plus: None,
            qual: None,
        }
    }

    /// # Usage
    /// Creates a new instance of a `[Record]` from a preinitialized
    /// `[Vec<u8>]` with the `id` and `seq` endpoints calculated. These
    /// endpoints are inclusive of the '\n' terminator and the data is
    /// expected to exclude the prefix '>' marker.
    /// ```
    /// let data = b">seq.0\nACGT\n".to_vec();
    /// let id = 6;
    /// let seq = 5;
    /// let fasta = fxread::Record::new_fasta(data, id, seq);
    /// assert_eq!(fasta.id(), b"seq.0");
    /// assert_eq!(fasta.seq(), b"ACGT");
    /// ```
    #[must_use]
    pub fn new_fasta(data: Vec<u8>, id: usize, seq: usize) -> Self {
        Self {
            data,
            id,
            seq,
            plus: None,
            qual: None,
        }
    }

    /// # Usage
    /// Creates a new instance of a `[Record]` from a preinitialized
    /// `[Vec<u8>]` with the `id`, `seq`, `plus`, and `qual` endpoints calculated.
    /// These endpoints are inclusive of the '\n' terminator and the data is
    /// expected to exclude the prefix '@' marker.
    /// ```
    /// let data = b"@seq.0\nACGT\n+\n1234\n".to_vec();
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
    #[must_use]
    pub fn new_fastq(data: Vec<u8>, id: usize, seq: usize, plus: usize, qual: usize) -> Self {
        Self {
            data,
            id,
            seq,
            plus: Some(plus),
            qual: Some(qual),
        }
    }

    /// Checks if `[Record]` is a fasta
    #[must_use]
    pub fn is_fasta(&self) -> bool {
        self.plus.is_none() & self.qual.is_none()
    }

    /// Checks if `[Record]` is a fastq
    #[must_use]
    pub fn is_fastq(&self) -> bool {
        self.plus.is_some() & self.qual.is_some()
    }

    /// Checks if `[Record]` has a valid header
    #[must_use]
    pub fn valid_header(&self) -> bool {
        if self.is_fasta() {
            self.data[0] == b'>'
        } else {
            self.data[0] == b'@'
        }
    }

    /// Checks if `[Record]` is empty
    #[must_use]
    pub fn empty(&self) -> bool {
        (self.id == 0) | (self.seq == 0)
    }

    /// Returns the Range of the ID
    #[must_use]
    pub fn id_range(&self) -> Range<usize> {
        1..self.id
    }

    /// Returns the Range of the sequence
    #[must_use]
    pub fn seq_range(&self) -> Range<usize> {
        1+self.id..self.id + self.seq
    }

    /// Returns the Range of the '+' region of a fastq
    #[must_use]
    pub fn plus_range(&self) -> Option<Range<usize>> {
        match self.plus {
            Some(plus) => Some(1+self.id + self.seq..self.id + self.seq + plus),
            None => None,
        }
    }

    /// Returns the Range of the quality score if it exists
    #[must_use]
    pub fn qual_range(&self) -> Option<Range<usize>> {
        let plus = match self.plus {
            Some(plus) => plus,
            None => return None,
        };
        match self.qual {
            Some(qual) => Some(1+self.id + self.seq + plus..self.id + self.seq + plus + qual),
            None => None,
        }
    }

    /// Returns a reference of the sequence ID
    #[must_use]
    pub fn id(&self) -> &[u8] {
        &self.data[self.id_range()]
    }

    /// Returns a reference of the sequence
    #[must_use]
    pub fn seq(&self) -> &[u8] {
        &self.data[self.seq_range()]
    }

    /// Returns a mutable reference of the sequence
    #[must_use]
    pub fn seq_mut(&mut self) -> &mut [u8] {
        let range = self.seq_range();
        &mut self.data[range]
    }

    /// Returns a reference of the '+' region of a fastq
    #[must_use]
    pub fn plus(&self) -> Option<&[u8]> {
        if let Some(range) = self.plus_range() {
            Some(&self.data[range])
        } else {
            None
        }
    }

    /// Returns a reference of the sequence
    #[must_use]
    pub fn qual(&self) -> Option<&[u8]> {
        if let Some(range) = self.qual_range() {
            Some(&self.data[range])
        } else {
            None
        }
    }

    /// Returns a mutable reference of the quality score if it exists
    #[must_use]
    pub fn qual_mut(&mut self) -> Option<&mut [u8]> {
        if let Some(range) = self.qual_range() {
            Some(&mut self.data[range])
        } else {
            None
        }
    }

    /// Returns a reference to the raw data underlying the record
    #[must_use]
    pub fn data(&self) -> &[u8] {
        &self.data
    }

    /// Validates that the record is not partially constructed
    /// or composed of unexpected nucleotides
    #[must_use]
    pub fn valid(&self) -> bool {
        if self.empty() {
            false
        } else {
            self.valid_sequence()
        }
    }

    /// Converts the sequence to uppercase
    #[must_use]
    pub fn seq_upper(&self) -> Vec<u8> {
        self.seq()
            .iter()
            .map(|c| if c & b' ' == 0 { *c } else { c ^ b' ' })
            .collect()
    }

    /// Reverse Complements the sequence
    #[must_use]
    pub fn seq_rev_comp(&self) -> Vec<u8> {
        self.seq()
            .iter()
            .rev()
            .map(|c| if c & 2 == 0 { c ^ 21 } else { c ^ 4 })
            .collect()
    }

    /// Converts all non-ACGTN nucleotides to N
    pub fn fix(&mut self) {
        self.seq_mut().iter_mut().for_each(|c| {
            if !matches!(
                c,
                b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'N' | b'n'
            ) {
                *c = b'N';
            }
        });
    }

    /// Converts the sequence to uppercase in place
    pub fn upper(&mut self) {
        self.seq_mut().iter_mut().for_each(|c| {
            if *c & b' ' != 0 {
                *c ^= b' ';
            }
        });
    }

    /// Reverse Complements the sequence in place
    /// Also reverses the quality scores if present
    pub fn rev_comp(&mut self) {
        // Reverse the sequence
        self.seq_mut().reverse();

        // Complement the sequence
        self.seq_mut().iter_mut().for_each(|c| {
            if *c & 2 == 0 {
                *c ^= 21;
            } else {
                *c ^= 4;
            }
        });

        // Reverse the quality scores if present
        if let Some(qual) = self.qual_mut() {
            qual.reverse();
        }
    }

    /// Validates whether sequence is composed
    /// of valid nucleotides
    fn valid_sequence(&self) -> bool {
        self.seq().iter().all(|b| {
            matches!(
                b,
                b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'N' | b'n' | b'U' | b'u'
            )
        })
    }

    /// Underlying record as str
    #[must_use]
    pub fn as_str_checked(&self) -> Result<&str, std::str::Utf8Error> {
        std::str::from_utf8(self.data())
    }

    /// ID as str
    #[must_use]
    pub fn id_str_checked(&self) -> Result<&str, std::str::Utf8Error> {
        std::str::from_utf8(self.id())
    }

    /// Sequence as str
    #[must_use]
    pub fn seq_str_checked(&self) -> Result<&str, std::str::Utf8Error> {
        std::str::from_utf8(self.seq())
    }

    /// Quality as str
    #[must_use]
    pub fn qual_str_checked(&self) -> Option<Result<&str, std::str::Utf8Error>> {
        if let Some(qual) = self.qual() {
            Some(std::str::from_utf8(qual))
        } else {
            None
        }
    }

    /// Underlying record as str unchecked (may panic if invalid utf8)
    #[must_use]
    pub fn as_str(&self) -> &str {
        self.as_str_checked().unwrap()
    }

    /// ID as str unchecked (may panic if invalid utf8)
    #[must_use]
    pub fn id_str(&self) -> &str {
        self.id_str_checked().unwrap()
    }

    /// Sequence as str unchecked (may panic if invalid utf8)
    #[must_use]
    pub fn seq_str(&self) -> &str {
        self.seq_str_checked().unwrap()
    }

    /// Quality as str unchecked (may panic if invalid utf8)
    #[must_use]
    pub fn qual_str(&self) -> Option<&str> {
        self.qual_str_checked().map(|qual| qual.unwrap())
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

impl Into<String> for Record {
    fn into(self) -> String {
        self.as_str().to_string()
    }
}

#[cfg(test)]
mod test {
    use super::Record;

    fn gen_valid_fasta() -> (Vec<u8>, usize, usize) {
        (b">seq.0\nACGT\n".to_vec(), 6, 5)
    }

    fn gen_invalid_fasta() -> (Vec<u8>, usize, usize) {
        (b">seq.0\nABCD\n".to_vec(), 6, 5)
    }

    fn gen_valid_fasta_lower() -> (Vec<u8>, usize, usize) {
        (b">seq.0\nacgt\n".to_vec(), 6, 5)
    }

    fn gen_invalid_fasta_lower() -> (Vec<u8>, usize, usize) {
        (b">seq.0\nabcd\n".to_vec(), 6, 5)
    }

    fn gen_valid_fasta_rev() -> (Vec<u8>, usize, usize) {
        (b">seq.0\nATCGGCTA\n".to_vec(), 6, 9)
    }

    fn gen_valid_fasta_rev_lower() -> (Vec<u8>, usize, usize) {
        (b">seq.0\natcggcta\n".to_vec(), 6, 9)
    }

    fn gen_valid_fastq() -> (Vec<u8>, usize, usize, usize, usize) {
        (b"@seq.0\nACGT\n+\n1234\n".to_vec(), 6, 5, 2, 5)
    }

    fn gen_invalid_fastq() -> (Vec<u8>, usize, usize, usize, usize) {
        (b"@seq.0\nABCD\n+\n1234\n".to_vec(), 6, 5, 2, 5)
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

    #[test]
    fn invalid_fix_fasta() {
        let (fasta, id, seq) = gen_invalid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(!record.valid());
        assert_eq!(record.seq(), b"ABCD");
        record.fix();
        assert!(record.valid());
        assert_eq!(record.seq(), b"ANCN");
    }

    #[test]
    fn invalid_fix_fastq() {
        let (fasta, id, seq, plus, qual) = gen_invalid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(!record.empty());
        assert!(!record.valid());
        assert_eq!(record.seq(), b"ABCD");
        record.fix();
        assert!(record.valid());
        assert_eq!(record.seq(), b"ANCN");
    }

    #[test]
    fn upper_inplace() {
        let (fasta, id, seq) = gen_valid_fasta_lower();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq(), b"acgt");
        record.upper();
        assert_eq!(record.seq(), b"ACGT");
    }

    #[test]
    fn upper_nochange_inplace() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq(), b"ACGT");
        record.upper();
        assert_eq!(record.seq(), b"ACGT");
    }

    #[test]
    fn reverse_complement_inplace() {
        let (fasta, id, seq) = gen_valid_fasta_rev();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq(), b"ATCGGCTA");
        record.rev_comp();
        assert_eq!(record.seq(), b"TAGCCGAT");
    }

    #[test]
    fn reverse_complement_fastq_inplace() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(!record.empty());
        assert!(record.valid());
        assert_eq!(record.seq(), b"ACGT");
        assert_eq!(record.qual().unwrap(), b"1234");
        record.rev_comp();
        assert_eq!(record.seq(), b"ACGT");
        assert_eq!(record.qual().unwrap(), b"4321");
    }

    #[test]
    fn fasta_str_methods() {
        let (fasta, id, seq) = gen_valid_fasta();
        let expected = String::from(">seq.0\nACGT\n");
        let record = Record::new_fasta(fasta, id, seq);
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.id_str_checked(), Ok("seq.0"));
        assert_eq!(record.seq_str(), "ACGT");
        assert_eq!(record.seq_str_checked(), Ok("ACGT"));
        let repr = record.as_str();
        assert_eq!(repr, &expected);
        let repr: String = record.into();
        assert_eq!(repr, expected);
    }

    #[test]
    fn fastq_str_methods() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let expected = String::from("@seq.0\nACGT\n+\n1234\n");
        let record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACGT");
        assert_eq!(record.qual_str(), Some("1234"));
        assert_eq!(record.id_str_checked(), Ok("seq.0"));
        assert_eq!(record.seq_str_checked(), Ok("ACGT"));
        assert_eq!(record.qual_str_checked(), Some(Ok("1234")));
        let repr = record.as_str();
        assert_eq!(repr, &expected);
        let repr: String = record.into();
        assert_eq!(repr, expected);
    }
}
