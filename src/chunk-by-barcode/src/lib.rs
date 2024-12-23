use std::fs::File;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{Read,BufReader,BufWriter};
use bio::io::fastq::{Reader,Writer,Record};
use itertools::{izip,Itertools};
use std::collections::HashMap;
use log::info;

/// Read a whitelist file into a sorted Vec
fn read_whitelist(file: &str) -> Vec<String> {
    let mut whitelist_file = File::open(file).unwrap();
    let mut whitelist: String = String::new();
    whitelist_file.read_to_string(&mut whitelist).unwrap();
    let mut whitelist: Vec<String> = whitelist.split("\n").map(|s| s.trim_end().to_string()).collect();
    // sort the whitelist, to ensure that regardless of whitelist sort order barcodes are always assigned to the same chunk for a given whitelist
    whitelist.sort();
    whitelist
}

/// Assign barcodes in a whitelist file to chunks
fn get_barcode_chunks(whitelist_file: &str, n_chunks: usize) -> HashMap<String,usize> {
    let mut whitelist = read_whitelist(whitelist_file);
    // add a '-' for reads w/o corrected barcode
    whitelist.push("-".to_string());
    let mut barcode_to_chunk: HashMap<String,usize> = HashMap::new();
    for (i, barcode) in whitelist.into_iter().enumerate() {
        let chunk = (i % n_chunks) + 1;
        barcode_to_chunk.insert(barcode, chunk);
    }
    barcode_to_chunk
}

/// Copy the description from one record to another
fn copy_description(from_read: &Record, to_read: &Record) -> Record {
    Record::with_attrs(
        to_read.id(),
        from_read.desc(),
        to_read.seq(),
        to_read.qual()
    )
}

/// Get the barcode for a record, based on the "CB" tag in the description
fn get_barcode(record: &Record) -> String {
    let description = record.desc().unwrap();
    let tags: Vec<&str> = description.split("\t").filter(|s| s.starts_with("CB:Z:")).collect();
    if tags.is_empty() {
        return "-".to_string();
    } else {
        assert!(tags.len() == 1);
        return tags[0].strip_prefix("CB:Z:").unwrap().to_string();
    }
}

struct Chunker {
    barcode_to_chunk: HashMap<String, usize>,
    writers: HashMap<usize, Writer<GzEncoder<File>>>,
}

impl Chunker {
    fn new(whitelist_file: &str, n_chunks: usize, prefix: &str) -> Self {
        let barcode_to_chunk = get_barcode_chunks(whitelist_file, n_chunks);
        let chunks: Vec<usize> = barcode_to_chunk.values().unique().cloned().collect();
        let mut writers: HashMap<usize, Writer<GzEncoder<File>>> = HashMap::new();
        for i in chunks {
            let filename = format!("{prefix}chunk_{i}.fastq.gz");
            let writer = Writer::from_bufwriter(BufWriter::new(GzEncoder::new(File::create(filename).unwrap(), Compression::fast())));
            writers.insert(i, writer);
        }

        Chunker {barcode_to_chunk, writers}
    }

    /// Write a fastq record to the correct chunk
    fn write_record(&mut self, record: &Record) {
        let barcode = get_barcode(record);
        let chunk = self.barcode_to_chunk.get(&barcode).expect(&format!("Could not find chunk for barcode: {barcode}"));
        let writer = self.writers.get_mut(chunk).unwrap();
        writer.write_record(record).unwrap();
    }

    /// Flush the writers.
    fn flush(&mut self) {
        for writer in self.writers.values_mut() {
            writer.flush().unwrap();
        }
    }
}


pub fn run(barcode_fastq: &str, read_1_fastq: &str, read_2_fastq: &str, whitelist_file: &str, n_chunks: usize) {
    let mut chunker_read_1 = Chunker::new(whitelist_file, n_chunks, "R1.");
    let mut chunker_read_2 = Chunker::new(whitelist_file, n_chunks, "R2.");

    let barcode_fastq_in = BufReader::new(MultiGzDecoder::new(File::open(barcode_fastq).unwrap()));
    let barcode_fastq_reader = Reader::from_bufread(barcode_fastq_in);

    let read_1_fastq_in = BufReader::new(MultiGzDecoder::new(File::open(read_1_fastq).unwrap()));
    let read_1_fastq_reader = Reader::from_bufread(read_1_fastq_in);

    let read_2_fastq_in = BufReader::new(MultiGzDecoder::new(File::open(read_2_fastq).unwrap()));
    let read_2_fastq_reader = Reader::from_bufread(read_2_fastq_in);

    let mut count: u64 = 0;

    for (barcode_result, read_1_result, read_2_result) in izip!(barcode_fastq_reader.records(), read_1_fastq_reader.records(), read_2_fastq_reader.records()) {

        count += 1;

        let barcode_read = barcode_result.unwrap();
        let read_1 = read_1_result.expect(&format!("Error at record {count}"));
        let read_2 = read_2_result.unwrap();

        // copy the description from barcode read to read 1 and read 2
        let new_read_1 = copy_description(&barcode_read, &read_1);
        let new_read_2 = copy_description(&barcode_read, &read_2);
        chunker_read_1.write_record(&new_read_1);
        chunker_read_2.write_record(&new_read_2);

        if count % 1_000_000 == 0 {
            info!("Processed {} read pairs so far", count);
        }
    }

    chunker_read_1.flush();
    chunker_read_2.flush();
}