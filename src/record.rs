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
        self.id.is_empty() & self.seq.is_empty()
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
}
