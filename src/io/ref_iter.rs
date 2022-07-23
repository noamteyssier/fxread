use std::io::Read;
use super::Parser;
use crate::record::{
    IdxRecord, 
    IdxRecordError,
    RefRecord
};
use anyhow::bail;

pub struct RecordRefIter<R: Read> {
    parser: Parser<R>,
    current: Option<IdxRecord>,
    current_length: Option<usize>
}
impl <R: Read> RecordRefIter <R> {

    pub fn new(p: Parser<R>) -> Self {
        Self {
            parser: p,
            current: None,
            current_length: None,
        }
    }

    pub fn next(&mut self) -> anyhow::Result<Option<RefRecord>> {
        match self.advance() {
            Ok(_) => Ok(self.get()),
            Err(e) => Err(e)
        }
    }

    pub fn get(&self) -> Option<RefRecord> {
        match self.current {
            None => None,
            Some(ref rec) => Some(rec.to_ref_record(self.parser.buffer.data()))
        }
    }

    fn step(is_fasta: bool, buffer: &[u8]) -> Result<IdxRecord, IdxRecordError> {
        match is_fasta {
            true => IdxRecord::fasta_from_buffer(buffer),
            false => IdxRecord::fastq_from_buffer(buffer),
        }

    }

    pub fn advance(&mut self) -> anyhow::Result<()> {
        let buffer = &mut self.parser.buffer;
        let mut reader = &mut self.parser.reader;

        if let Some(len) = self.current_length.take() {
            buffer.consume(len)
        }
        loop {
            match Self::step(self.parser.is_fasta, buffer.data()) {
                Ok(refrec) => {
                    self.current_length = Some(1 + refrec.end() - refrec.start());
                    self.current = Some(refrec);
                    break
                },

                Err(some_error) => match some_error {

                    // Buffer is empty (i.e. initialization or end)
                    IdxRecordError::EmptyBuffer => {
                        buffer.clean();
                        match buffer.read_from(&mut reader)? {
                            0 => {
                                self.current = None;
                                self.current_length = None;
                                break;
                            },
                            _ => { continue }
                        }
                    },

                    IdxRecordError::Truncated => {
                        buffer.clean();
                        if buffer.n_free() == 0 { bail!("Provided record is too long for current buffer") }
                        match buffer.read_from(&mut reader)? {
                            0 => { bail!("Unexpected truncation of file") },
                            _ => continue
                        }
                    }

                    IdxRecordError::MisplacedStart => { bail!("{:?}", some_error) },
                    IdxRecordError::MissingStart => { bail!("{:?}", some_error) },
                }
            }
        }

        Ok(())
    }
}

