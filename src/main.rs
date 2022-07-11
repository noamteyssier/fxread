use fxread::initialize_reader;

fn main() {
    let path = "sequences.fq";
    let reader = initialize_reader(path).expect("invalid path");
    reader
        .into_iter()
        .for_each(|x| println!("{:?}", x));
}
