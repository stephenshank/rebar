use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use bio::io::fastq;

use clap::{arg, Command};

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn main() {
    let matches = Command::new("rebar")
        .version("0.1.0")
        .author("Stephen D. Shank, Ph. D. <sshank@temple.edu>")
        .about("High performance barcode extractor")
        .arg(arg!(--fastq <VALUE>).required(true))
        .arg(arg!(--tsv <VALUE>).required(true))
        .get_matches();

    let fastq_file_path = matches.get_one::<String>("fastq").expect("required");
    println!("Input FASTQ: {:?}", fastq_file_path);
    println!(
        "Output TSV: {:?}",
        matches.get_one::<String>("tsv").expect("required")
    );
    let fastq_file = std::fs::File::open(fastq_file_path);
    let reader = fastq::Reader::new(fastq_file.unwrap());
    let barcode1 = b"AGTACGTACGAGTC";
    let barcode2 = b"GTACTCGCAGTAGTC";
    let scoring = Scoring {
        gap_open: -5,
        gap_extend: -1,
        match_fn: |a: u8, b: u8| if a == b { 1i32 } else { -3i32 },
        match_scores: Some((1, -3)),
        xclip_prefix: -10,
        xclip_suffix: MIN_SCORE,
        yclip_prefix: 0,
        yclip_suffix: 0,
    };

    let mut aligner = Aligner::with_capacity_and_scoring(1000, barcode2.len(), scoring);

    println!("Reading...");
    let mut read_vector: Vec<bio::io::fastq::Record> = Vec::new();
    let chunk_size = 5;
    for result in reader.records() {
        let record = result.unwrap();
        read_vector.push(record);
        let alignment1 = aligner.local(read_vector.last().unwrap().seq(), barcode1);
        println!("Aligning barcode 1...");
        println!(
            "{}",
            alignment1.pretty(read_vector.last().unwrap().seq(), barcode1)
        );
        println!("{},{}", alignment1.xstart, alignment1.xend);

        let alignment2 = aligner.local(read_vector.last().unwrap().seq(), barcode2);
        println!("Aligning barcode 2...");
        println!(
            "{}",
            alignment2.pretty(read_vector.last().unwrap().seq(), barcode2)
        );
        println!("{},{}", alignment2.xstart, alignment2.xend);
    }
}
