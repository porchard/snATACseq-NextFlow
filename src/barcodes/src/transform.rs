use std::fs::File;
use std::io::{Write,BufReader,BufWriter};
use std::collections::{HashSet,HashMap};
use std::fmt;
use log::info;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use bio::io::fastq;
use bio::alphabets::dna::revcomp;

use crate::whitelist;

/// Definition of the transformation needed to get a barcode out of a fastq sequence.
///
/// Has three elements
/// * the number of base pairs to remove from the start of the read
/// * the number of base pairs to remove from the end of the read
/// * whether to reverse complement the trimmed read
pub struct Transform {
    pub trim_from_start: usize,
    pub trim_from_end: usize,
    pub reverse_complement: bool,
}

impl fmt::Display for Transform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Transform(trim from start: {}; trim from end: {}; reverse_complement: {})", self.trim_from_start, self.trim_from_end, self.reverse_complement)
    }
}


/// Transform a fastq record.
///
/// # Example
/// ```
/// use bio::io::fastq::Record;
/// use barcodes::transform::Transform;
/// use barcodes::transform::transform_record;
/// let record = Record::with_attrs("read_name", Some("read_description"), b"CATGATGTTTTT", b"FAFFFFFFFFFF");
/// let transform_params = transform {trim_from_start: 1, trim_from_end: 5, reverse_complement: false};
/// let transformed_record = transform_record(&record, &transform_params);
/// assert_eq!(transformed_record, Record::with_attrs("read_name", Some("read_description"), b"ATGATG", b"AFFFFF"));
/// ```
pub fn transform_record(record: &fastq::Record, params: &Transform) -> fastq::Record {
    let trim_start = params.trim_from_start;
    let trim_end = record.seq().len() - params.trim_from_end;
    assert!(trim_start < trim_end);
    let trimmed_seq = &record.seq()[trim_start..trim_end];
    let trimmed_qual = &record.qual()[trim_start..trim_end];

    if params.reverse_complement {
    
        let rc_seq = revcomp(trimmed_seq);
        let r_qual: Vec<u8> = trimmed_qual.iter().rev().cloned().collect();
        let new_record = fastq::Record::with_attrs(
            record.id(), 
            record.desc(), 
            &rc_seq,
            &r_qual
        );

        new_record

    } else {
        
        let new_record = fastq::Record::with_attrs(
            record.id(), 
            record.desc(), 
            trimmed_seq,
            trimmed_qual
        );

        new_record
    }
}


fn infer_transform(fastq_filename: &str, whitelist_filename: &str, check_n_records: usize) -> Transform {

    let whitelist = whitelist::read_whitelist(whitelist_filename);
    let barcode_length = whitelist::barcode_length(&whitelist);
    let whitelist: HashSet<&[u8]> = whitelist.iter().map(|s| s.as_bytes()).collect();

    // read the first check_n_records records of the fastq file
    let mut transform_counts: Vec<(Transform, usize)> = Vec::new(); // this will store all the possible transforms that would yield a barcode of the correct length

    let fastq = BufReader::new(MultiGzDecoder::new(File::open(fastq_filename).unwrap()));
    let fastq_reader = fastq::Reader::from_bufread(fastq);
    let mut checked: usize = 0;

    for (i, result) in fastq_reader.records().enumerate() {

        let record = match result {
            Ok(x) => x,
            Err(x) => panic!("Error reading record {}: {}", i, x)
        };

        let read_length = record.seq().len();

        if i == 0 {
            // prepare the transforms
            for trim_start in 0..(read_length - barcode_length + 1) {
                let trim_end = trim_start + barcode_length;
                let trim_from_start = trim_start;
                let trim_from_end = read_length - trim_end;
                transform_counts.push((Transform{trim_from_start, trim_from_end, reverse_complement: true}, 0));
                transform_counts.push((Transform{trim_from_start, trim_from_end, reverse_complement: false}, 0));
            }
        }

        for (t, c) in transform_counts.iter_mut() {
            let transformed = transform_record(&record, t);
            if whitelist.contains(&transformed.seq()) {
                *c += 1;
            }
        }

        checked += 1;

        if i + 1 == check_n_records {
            break;
        }
    }

    for (t, c) in transform_counts.iter_mut() {
        info!("Potential transform {} matched whitelist {} of {} times", t, c, checked);
    }

    transform_counts.sort_by_key(|i| i.1);
    let (best_transform, count) = transform_counts.pop().unwrap();

    info!("Best transform: {} (matched whitelist {} of {} times)", best_transform, count, checked);

    best_transform
}


pub fn transform_fastq_file(input_fastq_filename: &str, whitelist_filename: &str, output_fastq_filename: &str, check_n_records: usize, output_counts_filename: &str) {

    let transform_params = infer_transform(input_fastq_filename, whitelist_filename, check_n_records);

    let whitelist = whitelist::read_whitelist(whitelist_filename);
    let whitelist: HashSet<&[u8]> = whitelist.iter().map(|s| s.as_bytes()).collect();

    let fastq_in = BufReader::new(MultiGzDecoder::new(File::open(input_fastq_filename).unwrap()));
    let fastq_reader = fastq::Reader::from_bufread(fastq_in);

    let fastq_out = BufWriter::new(GzEncoder::new(File::create(output_fastq_filename).unwrap(), Compression::fast()));
    let mut fastq_writer = fastq::Writer::from_bufwriter(fastq_out);

    let mut matched_whitelist: usize = 0;
    let mut total: usize = 0;
    let mut counts: HashMap<Vec<u8>,usize> = HashMap::new();

    for result in fastq_reader.records() {

        total += 1;

        let record = match result {
            Ok(x) => x,
            Err(x) => panic!("Error reading record {}: {}", total, x)
        };

        let transformed = transform_record(&record, &transform_params);
        fastq_writer.write_record(&transformed).unwrap();


        let count = counts.entry(transformed.seq().to_vec()).or_insert(0);
        *count += 1;
        
        if whitelist.contains(&transformed.seq()) {
            matched_whitelist += 1;
        }

        if total % 1_000_000 == 0 {
            info!("Processed {} reads so far ({} match whitelist)", total, matched_whitelist);
        }

    }

    fastq_writer.flush().unwrap();

    info!("Finished processing {} reads ({} match whitelist)", total, matched_whitelist);

    let mut counts_writer = BufWriter::new(File::create(output_counts_filename).unwrap());
    for (k, v) in counts.iter() {
        counts_writer.write(k).unwrap();
        counts_writer.write(b"\t").unwrap();
        counts_writer.write(v.to_string().as_bytes()).unwrap();
        counts_writer.write(b"\n").unwrap();
    }

    counts_writer.flush().unwrap();

}