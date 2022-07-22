use fxread::initialize_reader;

fn main() {
    let path = "mini.fq";
    let reader = initialize_reader(path).expect("invalid path");
    reader.count();
}
