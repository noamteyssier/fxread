use anyhow::{anyhow, Result};

#[derive(Debug)]
pub struct RefRecord <'a> {
    // (start, stop), but might include \r at the end
    head: usize,
    seq: usize,
    plus: Option<usize>,
    qual: Option<usize>,
    data: &'a [u8],
}
impl <'a> RefRecord <'a> {

    pub fn new(data: &'a [u8], head: usize, seq: usize) -> Self {

        Self {
            head, seq, data, plus: None, qual: None
        }

    }

    /// Returns a reference of the sequence ID
    pub fn id(&self) -> &[u8] {
        &self.data[1..self.head+1]
    }

    /// Returns a reference of the sequence 
    pub fn seq(&self) -> &[u8] {
        &self.data[self.head+2..self.head+self.seq+2]
    }

    pub fn str_id(&self) -> Result<&str> {
        match std::str::from_utf8(self.id()) {
            Ok(s) => Ok(s),
            Err(e) => Err(anyhow!("{:?}", e))
        }
    }

    pub fn str_seq(&self) -> Result<&str> {
        match std::str::from_utf8(self.seq()) {
            Ok(s) => Ok(s),
            Err(e) => Err(anyhow!("{:?}", e))
        }
    }

    pub fn idx_head(&self) -> usize {
        self.head
    }

    pub fn idx_seq(&self) -> usize {
        self.seq
    }

    pub fn data(&self) -> &[u8] {
        self.data
    }
}

