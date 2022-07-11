use std::io::BufRead;
use anyhow::Result;
use super::fastx::FastxRead;
use super::record::Record;

pub struct FastqReader {
    reader: Box<dyn BufRead>,
    buffer: String
}
impl FastqReader {
    pub fn new(reader: Box<dyn BufRead>) -> Self {
        Self { 
            reader,
            buffer: String::new()
        }
    }

    fn strip_header(&self, token: &str) -> String {
        token
            .trim_start_matches('@')
            .trim_end()
            .to_string()
    }

    fn strip_sequence(&self, token: &str) -> String {
        token
            .trim_end()
            .to_string()
    }
}

impl FastxRead for FastqReader {
    fn next_record(&mut self) -> Result<Option<Record>> {
        let mut record = Record::new();

        for index in 0..4 {
            if self.reader.read_line(&mut self.buffer)? == 0 { 
                break;
            }
            match index {
                0 => record.set_id(
                    self.strip_header(&self.buffer)
                ),
                1 => record.set_seq(
                    self.strip_sequence(&self.buffer)
                ),
                _ => ()
            };
            self.buffer.clear();
        }

        if record.empty() {
            Ok(None)
        }
        else {
            Ok(Some(record))
        }
    }

    fn print_records(&mut self) -> Result<()>{
        loop {
            match self.next_record()? {
                Some(r) => println!("{:?}", r),
                None => break
            };
        }
        Ok(())
    }
}

impl Iterator for FastqReader {

    type Item = Record;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(r) => r,
            Err(_) => panic!("Unexpected file end")
        }
    }

}

