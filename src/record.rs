use anyhow::{bail, Result};
use std::ops::{Range, RangeInclusive};

const DEFAULT_QUAL: u8 = b'F';

pub trait MyRange: Iterator<Item = i32> {
    fn start(&self) -> i32;
    fn end(&self) -> i32;
}
impl MyRange for Range<i32> {
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
}
impl MyRange for RangeInclusive<i32> {
    fn start(&self) -> i32 {
        *self.start()
    }
    fn end(&self) -> i32 {
        *self.end()
    }
}

/// An instance of a Fastx Record.
/// This is a two attribute object containing the sequence
/// ID and the Sequence.
#[derive(Debug)]
pub struct Record {
    data: Vec<u8>,
    /// The index of the ID size
    id: usize,
    /// The index of the sequence size
    seq: usize,
    /// The index of the plus size (if fastq)
    plus: Option<usize>,
    /// The index of the quality size (if fastq)
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

    /// # Usage
    ///
    /// Creates a new instance of a `[Record]` from its raw parts and calculates
    /// the endpoints. The `id` and `seq` are required and expected to be the raw
    /// data excluding the prefix '>' marker.
    ///
    /// ```
    /// let id = b"seq.0";
    /// let seq = b"ACGT";
    /// let record = fxread::Record::new_fasta_from_parts(id, seq).unwrap();
    /// assert_eq!(record.id(), b"seq.0");
    /// assert_eq!(record.seq(), b"ACGT");
    /// ```
    #[must_use]
    pub fn new_fasta_from_parts(id: &[u8], seq: &[u8]) -> Result<Self> {
        let mut data = Vec::with_capacity(id.len() + seq.len() + 2);
        if id.starts_with(b">") {
            bail!("ID cannot start with '>'");
        }
        if id.ends_with(b"\n") {
            bail!("ID cannot end with newline");
        }
        if seq.ends_with(b"\n") {
            bail!("Sequence cannot end with newline");
        }
        data.extend_from_slice(b">");
        data.extend_from_slice(id);
        data.extend_from_slice(b"\n");
        data.extend_from_slice(seq);
        data.extend_from_slice(b"\n");
        let id_idx = id.len() + 1;
        let seq_idx = seq.len() + 1;
        Ok(Self {
            data,
            id: id_idx,
            seq: seq_idx,
            plus: None,
            qual: None,
        })
    }

    /// # Usage
    ///
    /// Creates a new instance of a `[Record]` from its raw parts and calculates
    /// the endpoints. The `id` and `seq` and `qual` are required and expected to
    /// be the raw data excluding the prefix '>' or '@' markers. The `plus` is
    /// not required and is expected to be the raw data excluding the prefix '+'.
    ///
    /// ```
    /// let id = b"seq.0";
    /// let seq = b"ACGT";
    /// let qual = b"1234";
    /// let record = fxread::Record::new_fastq_from_parts(id, seq, qual).unwrap();
    /// assert_eq!(record.id(), b"seq.0");
    /// assert_eq!(record.seq(), b"ACGT");
    /// assert_eq!(record.qual().unwrap(), b"1234");
    /// ```
    #[must_use]
    pub fn new_fastq_from_parts(id: &[u8], seq: &[u8], qual: &[u8]) -> Result<Self> {
        let mut data = Vec::with_capacity(id.len() + seq.len() + qual.len() + 4);
        if id.starts_with(b"@") {
            bail!("ID cannot start with '@'");
        }
        if seq.len() != qual.len() {
            bail!("Sequence and Quality must be the same length");
        }
        if id.ends_with(b"\n") {
            bail!("ID cannot end with newline");
        }
        if seq.ends_with(b"\n") {
            bail!("Sequence cannot end with newline");
        }
        if qual.ends_with(b"\n") {
            bail!("Quality cannot end with newline");
        }
        data.extend_from_slice(b"@");
        data.extend_from_slice(id);
        data.extend_from_slice(b"\n");
        data.extend_from_slice(seq);
        data.extend_from_slice(b"\n+\n");
        data.extend_from_slice(qual);
        data.extend_from_slice(b"\n");
        let id_idx = id.len() + 1;
        let seq_idx = seq.len() + 1;
        let plus_idx = 2;
        let qual_idx = qual.len() + 1;
        Ok(Self {
            data,
            id: id_idx,
            seq: seq_idx,
            plus: Some(plus_idx),
            qual: Some(qual_idx),
        })
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
        1 + self.id..self.id + self.seq
    }

    /// Returns the Range of the '+' region of a fastq
    #[must_use]
    pub fn plus_range(&self) -> Option<Range<usize>> {
        match self.plus {
            Some(plus) => Some(1 + self.id + self.seq..self.id + self.seq + plus),
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
            Some(qual) => Some(1 + self.id + self.seq + plus..self.id + self.seq + plus + qual),
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

    /// Inserts nucleotides into the sequence at the specified index
    /// and the corresponding quality scores if present
    /// Returns an error if the index is greater than the sequence length
    pub fn insert_seq(&mut self, seq: &[u8], pos: usize) -> Result<()> {
        if pos >= self.seq {
            bail!("Cannot insert at a position greater than the sequence length");
        }
        let insert_pos = self.seq_range().start + pos;
        self.data
            .splice(insert_pos..insert_pos, seq.iter().cloned());
        self.seq += seq.len();
        self.insert_qual(seq.len(), pos)?;
        Ok(())
    }

    /// Inserts PHRED quality scores into the sequence at the specified index
    /// Returns an error if the index is greater than the qual length
    fn insert_qual(&mut self, insert_size: usize, pos: usize) -> Result<()> {
        if let Some(qual) = self.qual {
            if pos > qual {
                bail!("Cannot insert at a position greater than the quality length");
            }
            let insert_pos = self.qual_range().unwrap().start + pos;
            self.data
                .splice(insert_pos..insert_pos, vec![DEFAULT_QUAL; insert_size]);
            self.qual = Some(qual + insert_size);
        }
        Ok(())
    }

    /// Convenience function to insert nucleotides at the beginning of a sequence
    pub fn insert_seq_left(&mut self, seq: &[u8]) -> Result<()> {
        self.insert_seq(seq, 0)
    }

    /// Convenience function to insert nucleotides at the end of a sequence
    pub fn insert_seq_right(&mut self, seq: &[u8]) -> Result<()> {
        self.insert_seq(seq, self.seq - 1)
    }

    /// Trims the nucleotides from the left of the sequence
    /// and the corresponding quality scores if present
    /// Returns an error if the size is greater than the sequence length
    pub fn trim_left(&mut self, size: usize) -> Result<()> {
        if size > self.seq {
            bail!("Cannot trim more than the sequence length");
        }

        let seq_start = self.seq_range().start;
        let trim_end = seq_start + size;

        // Trim the sequence
        self.data.drain(seq_start..trim_end);

        // Update the seq index
        self.seq = self.seq.saturating_sub(size);

        // Fastq specific update (Trim quality scores)
        if self.is_fastq() {
            let qual_start = self.qual_range().unwrap().start;
            let trim_end = qual_start + size;
            self.data.drain(qual_start..trim_end);
            self.qual = Some(self.qual.unwrap().saturating_sub(size));
        }
        Ok(())
    }

    /// Trims the nucleotides from the right of the sequence
    /// and the corresponding quality scores if present
    /// Returns an error if the size is greater than the sequence length
    pub fn trim_right(&mut self, size: usize) -> Result<()> {
        if size > self.seq {
            bail!("Cannot trim more than the sequence length");
        }

        let seq_end = self.seq_range().end;
        let trim_start = seq_end - size;

        // Trim the sequence
        self.data.drain(trim_start..seq_end);

        // Update the seq index
        self.seq = self.seq.saturating_sub(size);

        // Fastq specific update (Trim quality scores)
        if self.is_fastq() {
            let qual_end = self.qual_range().unwrap().end;
            let trim_start = qual_end - size;
            self.data.drain(trim_start..qual_end);
            self.qual = Some(self.qual.unwrap().saturating_sub(size));
        }

        Ok(())
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

    #[test]
    fn fasta_trim_left_oversized() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(record.trim_left(10).is_err());
    }

    #[test]
    fn fasta_trim_right_oversized() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(record.trim_right(10).is_err());
    }

    #[test]
    fn fastq_trim_left_oversized() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(record.trim_left(10).is_err());
    }

    #[test]
    fn fastq_trim_right_oversized() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(record.trim_right(10).is_err());
    }

    #[test]
    fn fasta_trim_left() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        record.trim_left(2).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "GT");
        assert_eq!(record.as_str(), ">seq.0\nGT\n");
    }

    #[test]
    fn fastq_trim_left() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        record.trim_left(2).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "GT");
        assert_eq!(record.qual_str(), Some("34"));
        assert_eq!(record.as_str(), "@seq.0\nGT\n+\n34\n");
    }

    #[test]
    fn fasta_trim_right() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        record.trim_right(2).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "AC");
        assert_eq!(record.as_str(), ">seq.0\nAC\n");
    }

    #[test]
    fn fastq_trim_right() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        record.trim_right(2).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "AC");
        assert_eq!(record.qual_str(), Some("12"));
        assert_eq!(record.as_str(), "@seq.0\nAC\n+\n12\n");
    }

    #[test]
    fn fasta_insert_left() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        record.insert_seq_left(b"TT").unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "TTACGT");
        assert_eq!(record.as_str(), ">seq.0\nTTACGT\n");
    }

    #[test]
    fn fastq_insert_left() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        record.insert_seq_left(b"TT").unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "TTACGT");
        assert_eq!(record.qual_str(), Some("FF1234"));
        assert_eq!(record.as_str(), "@seq.0\nTTACGT\n+\nFF1234\n");
    }

    #[test]
    fn fasta_insert_right() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        record.insert_seq_right(b"TT").unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACGTTT");
        assert_eq!(record.as_str(), ">seq.0\nACGTTT\n");
    }

    #[test]
    fn fastq_insert_right() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        record.insert_seq_right(b"TT").unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACGTTT");
        assert_eq!(record.qual_str(), Some("1234FF"));
        assert_eq!(record.as_str(), "@seq.0\nACGTTT\n+\n1234FF\n");
    }

    #[test]
    fn fasta_insert_middle() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        record.insert_seq(b"TT", 2).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACTTGT");
        assert_eq!(record.as_str(), ">seq.0\nACTTGT\n");
    }

    #[test]
    fn fastq_insert_middle() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        record.insert_seq(b"TT", 2).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACTTGT");
        assert_eq!(record.qual_str(), Some("12FF34"));
        assert_eq!(record.as_str(), "@seq.0\nACTTGT\n+\n12FF34\n");
    }

    #[test]
    fn fasta_insert_custom_oversized() {
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(record.insert_seq(b"TT", 10).is_err());

        // seq size is error because it includes the newline
        let (fasta, id, seq) = gen_valid_fasta();
        let mut record = Record::new_fasta(fasta, id, seq);
        assert!(record.insert_seq(b"TT", seq).is_err());
    }

    #[test]
    fn fastq_insert_custom_oversized() {
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(record.insert_seq(b"TT", 10).is_err());

        // seq size is error because it includes the newline
        let (fasta, id, seq, plus, qual) = gen_valid_fastq();
        let mut record = Record::new_fastq(fasta, id, seq, plus, qual);
        assert!(record.insert_seq(b"TT", seq).is_err());
    }

    #[test]
    fn init_fastq_from_parts() {
        let id = b"seq.0";
        let seq = b"ACGT";
        let qual = b"1234";
        let record = Record::new_fastq_from_parts(id, seq, qual).unwrap();
        let expected_plus = vec![b'+'];
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACGT");
        assert_eq!(record.plus(), Some(expected_plus.as_slice()));
        assert_eq!(record.qual_str(), Some("1234"));
    }

    #[test]
    fn init_fasta_from_parts() {
        let id = b"seq.0";
        let seq = b"ACGT";
        let record = Record::new_fasta_from_parts(id, seq).unwrap();
        assert_eq!(record.id_str(), "seq.0");
        assert_eq!(record.seq_str(), "ACGT");
        assert!(record.plus().is_none());
        assert!(record.qual_str().is_none());
    }

    #[test]
    fn init_fastq_from_parts_invalid_id_a() {
        let id = b"seq.0\n";
        let seq = b"ACGT";
        let qual = b"1234";
        assert!(Record::new_fastq_from_parts(id, seq, qual).is_err());
    }

    #[test]
    fn init_fastq_from_parts_invalid_id_b() {
        let id = b"@seq.0";
        let seq = b"ACGT";
        let qual = b"1234";
        assert!(Record::new_fastq_from_parts(id, seq, qual).is_err());
    }

    #[test]
    fn init_fastq_from_parts_invalid_seq_a() {
        let id = b"seq.0";
        let seq = b"ACGT\n";
        let qual = b"1234";
        assert!(Record::new_fastq_from_parts(id, seq, qual).is_err());
    }

    #[test]
    fn init_fastq_from_parts_invalid_seq_qual_size_mismatch() {
        let id = b"seq.0";
        let seq = b"ACGTA";
        let qual = b"1234";
        assert!(Record::new_fastq_from_parts(id, seq, qual).is_err());
    }

    #[test]
    fn init_fastq_from_parts_invalid_qual_a() {
        let id = b"seq.0";
        let seq = b"ACGT";
        let qual = b"1234\n";
        assert!(Record::new_fastq_from_parts(id, seq, qual).is_err());
    }

    #[test]
    fn init_fasta_from_parts_invalid_id_a() {
        let id = b"seq.0\n";
        let seq = b"ACGT";
        assert!(Record::new_fasta_from_parts(id, seq).is_err());
    }

    #[test]
    fn init_fasta_from_parts_invalid_id_b() {
        let id = b">seq.0";
        let seq = b"ACGT";
        assert!(Record::new_fasta_from_parts(id, seq).is_err());
    }

    #[test]
    fn init_fasta_from_parts_invalid_seq_a() {
        let id = b"seq.0";
        let seq = b"ACGT\n";
        assert!(Record::new_fasta_from_parts(id, seq).is_err());
    }
}
