use memchr::memchr;
use super::{
    RefRecord,
    Record
};

#[derive(Debug)]
pub enum IdxRecordError {
    EmptyBuffer,
    MisplacedStart,
    MissingStart,
    Truncated
}

#[derive(Debug)]
pub struct IdxRecord {
    head: usize,
    seq: usize,
    plus: Option<usize>,
    qual: Option<usize>,
    start: usize,
    end: usize
}

impl IdxRecord {

    fn new_fasta(head: usize, seq: usize, start: usize, end: usize) -> Self {
        Self {
            head, seq, start, end,
            plus: None, qual: None
        }
    }

    fn new_fastq(head: usize, seq: usize, plus: usize, qual: usize, start: usize, end: usize) -> Self {
        Self {
            head, seq, start, end,
            plus: Some(plus), qual: Some(qual)
        }
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }
    
    fn has_marker(marker: u8, buffer: &[u8]) -> Result<(), IdxRecordError> {
        match memchr(marker, buffer) {
            Some(index) => match index {
                0 => Ok(()),
                _ => Err(IdxRecordError::MisplacedStart)
            },
            None => Err(IdxRecordError::MissingStart)

        }
    }

    fn capture_newline(buffer: &[u8]) -> Result<usize, IdxRecordError> {
        match memchr(b'\n', buffer) {
            Some(index) => Ok(index),
            None => Err(IdxRecordError::Truncated)
        }
    }

    pub fn fasta_from_buffer(buffer: &[u8]) -> Result<Self, IdxRecordError> {
        if buffer.len() == 0 { return Err(IdxRecordError::EmptyBuffer) }

        Self::has_marker(b'>', buffer)?;
        let head = Self::capture_newline(&buffer[1..])?;
        let seq = Self::capture_newline(&buffer[head+2..])?;
        Ok( 
            Self::new_fasta(head, seq, 0, head+seq+2)    
        )
    }

    pub fn fastq_from_buffer(buffer: &[u8]) -> Result<Self, IdxRecordError> {
        if buffer.len() == 0 { return Err(IdxRecordError::EmptyBuffer) }

        Self::has_marker(b'@', buffer)?;
        let head = Self::capture_newline(&buffer[1..])?;
        let seq = Self::capture_newline(&buffer[head+2..])?;
        let plus = Self::capture_newline(&buffer[head+seq+3..])?;
        let qual = Self::capture_newline(&buffer[head+seq+plus+4..])?;

        Ok(
            Self::new_fastq(head, seq, plus, qual, 0, head+seq+plus+qual+4)
        )
    }

    pub fn to_ref_record<'a>(&self, buffer: &'a [u8]) -> RefRecord<'a> {
        let data = &buffer[self.start..self.end];
        RefRecord::new(data, self.head, self.seq)
    }

    pub fn to_record(&self, buffer: &[u8]) -> Record {
        let data = &buffer[self.start..self.end];
        Record::new(data, self.head, self.seq)
    }
}
