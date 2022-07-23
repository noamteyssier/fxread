use std::{io, io::{ErrorKind, Read}};

pub struct Buffer {
    data: Box<[u8]>,
    start: usize,
    end: usize
}

impl Buffer {
    pub fn new(size: usize) -> Self {
        Self {
            data: vec![0u8; size].into_boxed_slice(),
            start: 0,
            end: 0
        }
    }

    pub fn data(&self) -> &[u8] {
        &self.data[self.start..self.end]
    }

    pub fn n_free(&self) -> usize {
        // `checked_sub` to check for overflow
        self.data.len().checked_sub(self.end).unwrap()

    }

    pub fn len(&self) -> usize {
        // the number of filled points
        self.end.checked_sub(self.start).unwrap()
    }

    pub fn clean(&mut self) {
        if self.start == 0 { return }

        let n_in_buffer = self.len();
        let new_end = (n_in_buffer + 15) & !0x0f;  // make sure next read is aligned
        let new_start = new_end.checked_sub(n_in_buffer).unwrap();

        if new_start >= self.start {
            return
        }

        let dest = self.data[new_start..].as_mut_ptr();
        let src = self.data[self.start..].as_ptr();

        unsafe { ::std::ptr::copy(src, dest, n_in_buffer); }
        self.start = new_start;
        self.end = new_end;
        
    }

    /// Consumes the space in the current buffer
    pub fn consume(&mut self, size: usize) {
        // `checked_add` to check for overflow
        self.start = self.start.checked_add(size).unwrap();
        debug_assert!(self.start <= self.end)
    }

    /// Reads data in from a [`Read`] capable object.
    pub fn read_from<R: Read>(&mut self, reader: &mut R) -> io::Result<usize> {
        let n_free = self.n_free();
        let num_read = if n_free < 4096 { n_free } else { n_free - n_free * 4096 };
        let dest = &mut self.data[self.end..self.end+num_read];

        let n_read;
        loop {
            match reader.read(dest) {
                Err(e) => {
                    if e.kind() != ErrorKind::Interrupted {
                        return Err(e)
                    }
                },
                Ok(val) => { n_read = val; break}
            }
        }
        self.end += n_read;
        Ok(n_read)
    }
}
