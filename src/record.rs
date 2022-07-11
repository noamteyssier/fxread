
#[derive(Debug)]
pub struct Record {
    id: String,
    seq: String
}
impl Record {
    pub fn new() -> Self {
        Self {
            id: String::new(),
            seq: String::new()
        }
    }

    pub fn empty(&self) -> bool {
        self.id.is_empty() & self.seq.is_empty()
    }

    pub fn set_id(&mut self, token: String) {
        self.id = token;
    }

    pub fn set_seq(&mut self, token: String) {
        self.seq = token;
    }
}
