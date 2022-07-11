use super::Record;
use anyhow::Result;

pub trait FastxRead : Iterator {
    fn next_record(&mut self) -> Result<Option<Record>>;
    fn print_records(&mut self) -> Result<()>;
}
