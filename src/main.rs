use clap::{arg, Command};

fn main() {
    let matches = Command::new("rustybarcoder")
        .version("0.1.0")
        .author("Stephen D. Shank, Ph. D. <sshank@temple.edu>")
        .about("High performance barcode extractor")
        .arg(arg!(--fastq <VALUE>).required(true))
        .arg(arg!(--tsv <VALUE>).required(true))
        .get_matches();

    println!(
        "Input FASTQ: {:?}",
        matches.get_one::<String>("fastq").expect("required")
    );
    println!(
        "Output TSV: {:?}",
        matches.get_one::<String>("tsv").expect("required")
    );
}
