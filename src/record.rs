use anyhow::{anyhow, Result};

/// An instance of a Fastx Record.
/// This is a two attribute object containing the sequence
/// ID and the Sequence.
#[derive(Debug)]
pub struct Record {
    id: String,
    seq: String
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
            id: String::new(),
            seq: String::new()
        }
    }

    /// Checks if `[Record]` is empty
    pub fn empty(&self) -> bool {
        self.id.is_empty() | self.seq.is_empty()
    }

    /// # Usage
    /// Sets the ID of the record
    /// ```
    /// let mut record = fxread::Record::new();
    /// record.set_id("some_id".to_string());
    /// assert_eq!(record.id(), "some_id");
    /// ```
    pub fn set_id(&mut self, token: String) {
        self.id = token;
    }

    /// # Usage
    /// Sets the Sequence of the record
    /// ```
    /// let mut record = fxread::Record::new();
    /// record.set_seq("ACGT".to_string());
    /// assert_eq!(record.seq(), "ACGT");
    /// ```
    pub fn set_seq(&mut self, token: String) {
        self.seq = token;
    }

    /// Returns a reference of the sequence ID
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns a reference of the sequence 
    pub fn seq(&self) -> &str {
        &self.seq
    }

    /// Validates that the record is not partially constructed
    /// or composed of unexpected nucleotides
    pub fn valid(&self) -> bool {
        !self.empty() & self.valid_sequence()
    }

    /// Converts the sequence to uppercase
    pub fn seq_upper(&self) -> String {
        self.seq.to_ascii_uppercase()
    }

    /// Reverse Complements the sequence
    /// This will also convert the sequence to uppercase
    pub fn seq_rev_comp(&self) -> Result<String> {
        self.seq
            .chars()
            .rev()
            .map(|c| match c {
                'A'|'a' => Ok('T'),
                'C'|'c' => Ok('G'),
                'G'|'g' => Ok('C'),
                'T'|'t' => Ok('A'),
                'N'|'n' => Ok('N'),
                _ => Err(anyhow!("Unexpected nucleotide found: {}", c))
            }).collect()
    }

    /// Validates whether sequence is composed
    /// of valid nucleotides
    fn valid_sequence(&self) -> bool {
        self.seq.chars().all(|c| match c {
            'A'|'a'|'C'|'c'|'G'|'g'|'T'|'t'|'N'|'n'|'U'|'u' => true,
            _ => false
        })
    }

}

#[cfg(test)]
mod test {
    use super::Record;

    #[test]
    fn create() {
        let record = Record::new();
        assert!(record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn create_partial_id() {
        let mut record = Record::new();
        record.set_id(String::from("some_id"));
        assert!(record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn create_partial_seq() {
        let mut record = Record::new();
        record.set_seq(String::from("ACGT"));
        assert!(record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn valid() {
        let mut record = Record::new();
        record.set_id(String::from("some_id"));
        record.set_seq(String::from("ACGT"));
        assert!(!record.empty());
        assert!(record.valid());
    }

    #[test]
    fn invalid() {
        let mut record = Record::new();
        record.set_id(String::from("some_id"));
        record.set_seq(String::from("BCGT"));
        assert!(!record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn valid_lowercase() {
        let mut record = Record::new();
        record.set_id(String::from("some_id"));
        record.set_seq(String::from("acgt"));
        assert!(!record.empty());
        assert!(record.valid());
    }

    #[test]
    fn invalid_lowercase() {
        let mut record = Record::new();
        record.set_id(String::from("some_id"));
        record.set_seq(String::from("bcgt"));
        assert!(!record.empty());
        assert!(!record.valid());
    }

    #[test]
    fn upper_conversion() {
        let mut record = Record::new();
        record.set_seq(String::from("acgt"));
        assert_eq!(record.seq_upper(), String::from("ACGT"));
    }

    #[test]
    fn reverse_complement() {
        let mut record = Record::new();
        record.set_seq(String::from("ACGTA"));
        assert_eq!(record.seq_rev_comp().unwrap(), "TACGT");
    }

    #[test]
    fn lower_reverse_complement() {
        let mut record = Record::new();
        record.set_seq(String::from("acgta"));
        assert_eq!(record.seq_rev_comp().unwrap(), "TACGT");
    }

    #[test]
    #[should_panic]
    fn invalid_reverse_complement() {
        let mut record = Record::new();
        record.set_seq(String::from("ACGTAB"));
        record.seq_rev_comp().unwrap();
    }
}
