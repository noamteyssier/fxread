#![feature(test)]

use std::{io::BufReader, fs::File};
use flate2::read::MultiGzDecoder;
use anyhow::Result;

extern crate test;
use test::Bencher;

mod record;
mod fastx;
mod fastqread_generic;
mod fastqread_dyn;

use fastx::FastxRead;
use record::Record;



fn initialize_generic_reader(path: &str) -> Result<Box<dyn FastxRead<Item = Record>>> {
    match path.ends_with(".gz") {
        true => {
            let file = File::open(path)?;
            let gzip = MultiGzDecoder::new(file);
            let buffer = BufReader::new(gzip);
            let reader = fastqread_generic::FastqReader::new(buffer);
            Ok(Box::new(reader))
        },
        false => {
            let file = File::open(path)?;
            let buffer = BufReader::new(file);
            let reader = fastqread_generic::FastqReader::new(buffer);
            Ok(Box::new(reader))
        }
    }
}

fn initialize_dyn_reader(path: &str) -> Result<fastqread_dyn::FastqReader> {
    match path.ends_with(".gz") {
        true => {
            let file = File::open(path)?;
            let gzip = MultiGzDecoder::new(file);
            let buffer = BufReader::new(gzip);
            let reader = fastqread_dyn::FastqReader::new(Box::new(buffer));
            Ok(reader)
        },
        false => {
            let file = File::open(path)?;
            let buffer = BufReader::new(file);
            let reader = fastqread_dyn::FastqReader::new(Box::new(buffer));
            Ok(reader)
        }
    }
}

#[bench]
fn benchmark_generic(b: &mut Bencher){
    let path = "sequences.fq.gz";

    b.iter(
        || {
            (0..100).for_each(|_| {
                let reader = initialize_generic_reader(path).unwrap();
                reader
                    .into_iter()
                    .for_each(|_| {});
            })
        });
}

#[bench]
fn benchmark_dynamic(b: &mut Bencher){
    let path = "sequences.fq.gz";
    b.iter(
        || {
            (0..100).for_each(|_| {
                let reader = initialize_dyn_reader(path).unwrap();
                reader
                    .into_iter()
                    .for_each(|_| {});
            })
        });
}

fn main() -> Result<()> {
    Ok(())
}
