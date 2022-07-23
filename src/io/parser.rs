use std::io::Read;
use super::Buffer;
use super::RecordRefIter;

const BUFSIZE: usize = 1024 * 1;

pub struct Parser<R: Read> {
    pub reader: R,
    pub buffer: Buffer,
    pub is_fasta: bool 
}
impl <R: Read> Parser <R> {
    pub fn new_fasta(reader: R) -> Self {
        Self {
            reader, 
            buffer: Buffer::new(BUFSIZE), 
            is_fasta: true
        }
    }

    pub fn new_fastq(reader: R) -> Self {
        Self {
            reader, 
            buffer: Buffer::new(BUFSIZE), 
            is_fasta: false
        }
    }

    pub fn ref_iter(self) -> RecordRefIter<R> {
        RecordRefIter::new(self)
    }
}
