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


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_fasta() {
        let path = "example/sequences.fa";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn test_gzfasta() {
        let path = "example/sequences.fa.gz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn test_fastq() {
        let path = "example/sequences.fq";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }

    #[test]
    fn test_gzfastq() {
        let path = "example/sequences.fq.gz";
        let reader = initialize_reader(path).expect("invalid path");
        let num_records = reader
            .into_iter()
            .map(|x| assert!(!x.empty()))
            .count();
        assert_eq!(num_records, 10);
    }
}
