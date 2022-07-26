# fxread
A barebones fastx reader for rust.

This crate attempts to be a faster and more lightweight alternative to `bio-rs` and provides a standardized interface to working with fasta and fastq formats.
The goal of this crate is to be fast and flexible - it is about twice as fast as `bio-rs` on average but about half as fast than `fastq` for standard fastx files. 
The difference between the different crates is reduced heavily though once gzip files are included. 

The speed up can be attributed to reducing the total number of vectors allocated for each record - but the limitation compared to `fastq` is that each record has ownership over its data and is allocated once.
This creates extra overhead, but is very convenient as you can treat the reader directly as an iterator.

Some benefits of this interface is that each `FastaReader` and `FastqReader` share the `FastxReader` trait and act as iterators over `Record`s.
An example of the interface:
```rust
use fxread::initialize_reader;
let path = "example/sequences.fq";
let reader = initialize_reader(path).unwrap();
assert_eq!(reader.count(), 10);
```

Check out the [API Documentation](https://docs.rs/fxread) for usage
