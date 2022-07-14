pub mod record;
pub mod fastx;
pub mod fastq;
pub mod fasta;
pub mod utils;

pub use fastx::FastxRead;
pub use fasta::FastaReader;
pub use fastq::FastqReader;
pub use record::Record;
pub use utils::initialize_reader;
