//! # Summary
//! This is a package to handle reading of fastx records.
//!
//! fastx records are defined from fasta or fastq formats but a fastx
//! record is one with an id and sequence regardless of its input format.
//!
//! This package defines the [`FastxRead`] trait, which is implemented by both
//! [`FastaReader`] and [`FastqReader`]. This trait allows these types to be
//! treated as an [`Iterator`] which returns a fastx [`Record`] on each `next()`.
//!
//! This package also handles the creation of a [`FastxRead`] capable reader from
//! the naming of an input file using dynamic dispatch. This utility is found in
//! [`initialize_reader`]. Please see usage for example usages or each items
//! unit tests.
//!
//! # Usage
//! ## Reading from a file directly
//! This is a very common usecase, where you have some file
//! you'd like to read records from directly. Check out this
//! function directly for more examples.
//! ```
//! use fxread::initialize_reader;
//! let path = "example/sequences.fa";
//! let reader = initialize_reader(path).unwrap();
//! reader
//!     .for_each(|record| println!("{:?}", record));
//! ```

/// Module for fastx records.
pub mod record;

/// Module for the fastxread trait
pub mod fastx;

/// Module for a fastq reader
pub mod fastq;

/// Module for a fasta reader
pub mod fasta;

/// Module for utility functions associated with creating
/// the correct fastx reader.
pub mod utils;

pub use fasta::FastaReader;
pub use fastq::FastqReader;
pub use fastx::FastxRead;
pub use record::Record;
pub use utils::{initialize_reader, initialize_stdin_reader};
