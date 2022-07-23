use super::RefRecord;

pub struct Record {
    head: usize,
    seq: usize,
    data: Vec<u8>
}
impl Record {
    pub fn new(data: &[u8], head: usize, seq: usize) -> Self {
        Self {
            data: data.to_owned(),
            head,
            seq
        }
    }
    pub fn from_reference(refrec: RefRecord) -> Self {
        Self {
            head: refrec.idx_head(),
            seq: refrec.idx_seq(),
            data: refrec.data().to_owned()
        }
    }
}

