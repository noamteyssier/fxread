use fxread::initialize_reader;

fn main() {
    let path = "example/sequences.fa";
    let reader = initialize_reader(path).expect("invalid path");
    println!("{}", reader.count());
}
