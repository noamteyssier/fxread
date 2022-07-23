use std::io::Read;
use std::fs::File;
use fxread::io::Parser;

fn reference_iter<R: Read>(reader: R) -> anyhow::Result<()> {
    let parser = Parser::new_fastq(reader);
    let mut iter = parser.ref_iter();
    let mut count = 0;
    while let Some(_rec) = iter.next()? {
        count += 1;
    }
    println!("{}", count);
    Ok(())
}


fn main() {
    //let path = "example/sequences.fa";
    let path = "mini.fq";
    let file = File::open(path).unwrap();
    reference_iter(file).unwrap();
    //let reader = initialize_reader(path).expect("invalid path");
    //reader.count();
}
