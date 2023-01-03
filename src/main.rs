use bio::alignment::pairwise::*;
use bio::io::fastq;
use rayon::prelude::*;
use std::fs::File;
use std::fs::OpenOptions;

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
    let tsv_file_path = matches.get_one::<String>("tsv").expect("required");
    println!("Input FASTQ: {:?}", fastq_file_path);
    println!("Output TSV: {:?}", tsv_file_path);
    let fastq_file = std::fs::File::open(fastq_file_path);
    let reader = fastq::Reader::new(fastq_file.unwrap());
    let barcode1 = b"AGTACGTACGAGTC";
    let barcode2 = b"GTACTCGCAGTAGTC";

    println!("Reading...");

    let tsv_file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(tsv_file_path)
        .unwrap();
    let mut read_vector: Vec<bio::io::fastq::Record> = Vec::new();
    let mut wtr = csv::Writer::from_writer(tsv_file);
    wtr.write_record(&["b1start", "b2start"]);
    let chunk_size = 10000;
    for (read_index, result) in reader.records().enumerate() {
        let record = result.unwrap();
        read_vector.push(record);
        if read_index.rem_euclid(chunk_size) == 0 {
            let result_vector: Vec<(usize, usize)> = read_vector
                .par_iter()
                .map(|result| {
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

                    let mut aligner =
                        Aligner::with_capacity_and_scoring(1000, barcode2.len(), scoring);
                    let alignment1 = aligner.local(result.seq(), barcode1);
                    let alignment2 = aligner.local(result.seq(), barcode2);
                    (alignment1.xstart, alignment2.xstart)
                })
                .collect();
            for result in result_vector {
                let (b1_start, b2_start) = result;
                wtr.write_record([b1_start.to_string(), b2_start.to_string()]);
            }
            read_vector.clear();
        }
    }
}
