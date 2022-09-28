use super::Record;
use anyhow::Result;

/// A trait for Fasta and Fastq readers
pub trait FastxRead: Iterator {
    /// Returns the next fastx [`Record`] in the iterator.
    fn next_record(&mut self) -> Result<Option<Record>>;
}
