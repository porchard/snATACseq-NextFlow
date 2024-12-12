use std::fs::File;
use std::collections::{HashSet,HashMap};
use std::cmp;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{Read,BufReader,BufWriter};
use bio::io::fastq;
use itertools::izip;
use log::info;
use crate::trie::Trie;
use crate::whitelist;


struct BarcodeCorrector {
    whitelist: HashSet<String>,
    counts: HashMap<String,usize>,
    whitelist_trie: Trie,
}

impl BarcodeCorrector {
    pub fn from(whitelist_barcodes: &HashSet<String>, barcode_counts: &HashMap<String,usize>) -> Self {
        let pseudocount: usize = 1;

        let whitelist: HashSet<String> = whitelist_barcodes.iter().cloned().collect();
        let mut counts: HashMap<String, usize> = HashMap::new();
        let mut whitelist_trie = Trie::new();
        
        for i in whitelist.iter() {
            whitelist_trie.add_word(i.as_bytes());
            counts.insert(i.clone(), barcode_counts.get(i).copied().unwrap_or(0) + pseudocount);
        }
        
        BarcodeCorrector { whitelist, counts, whitelist_trie }
    }


    pub fn correct_barcode(&self, barcode: &[u8], phred: &[u8], max_distance: usize) -> Option<String> {

        let uncorrected = String::from_utf8(barcode.to_vec()).unwrap();

        if self.whitelist.contains(&uncorrected) {
            return Some(uncorrected);
        }

        let mut similar: Vec<String> = self.whitelist_trie.get_words_within_hamming_distance(barcode, max_distance).into_iter().map(|(i, _j)| i).collect();

        if similar.is_empty() {
            return None;
        } else if similar.len() == 1 {
            return similar.pop();
        } else {
            let probability_of_errors: Vec<f64> = similar.iter().map(|s| probability_of_incorrect_base_calls(barcode, s.as_bytes(), phred)).collect();
            let probability_of_errors_times_count: Vec<f64> = probability_of_errors.into_iter().enumerate().map(|(i, p)| p*(*self.counts.get(&similar[i]).unwrap() as f64)).collect();
            let norm_factor: f64 = probability_of_errors_times_count.iter().sum();
            let posteriors: Vec<f64> = probability_of_errors_times_count.iter().map(|i| i / norm_factor).collect();
    
            for (correction, p) in izip!(similar, posteriors) {
                if p >= 0.975 {
                    return Some(correction);
                }
            }
        }
    
        None
    
    }
}



// To match CellRanger corrections, max_allowed_quality should be 66
fn probability_of_incorrect_base_call(quality_score: &u8, max_quality_score: &u8) -> f64 {
    let q = (cmp::min(*quality_score, *max_quality_score) as f64) - 33.0;
    let power_base: f64 = 10.0;
    power_base.powf(-1.0 * q / 10.0)
}

fn probability_of_incorrect_base_calls(uncorrected: &[u8], corrected: &[u8], phred: &[u8]) -> f64 {

    assert_eq!(uncorrected.len(), corrected.len());
    
    let mut l: f64 = 1.0;
    
    for (u, c, p) in izip!(uncorrected, corrected, phred) {
        if u != c {
            l *= probability_of_incorrect_base_call(p, &66);
        }
    }

    l
}


fn read_barcode_counts(file: &str) -> HashMap<String, usize> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut counts_file = File::open(file).unwrap();
    let mut counts_string = String::new();
    counts_file.read_to_string(&mut counts_string).unwrap();
    counts_string = counts_string.trim().to_string();
    for i in counts_string.split("\n") {
        let barcode_and_count: Vec<&str> = i.split("\t").collect();
        let barcode = barcode_and_count[0].to_string();
        let count = barcode_and_count[1].parse::<usize>().unwrap();
        let e = counts.entry(barcode).or_insert(0);
        *e += count;
    }

    counts
}


pub fn correct_barcodes_in_fastq(input_fastq_filename: &str, whitelist_filename: &str, counts_filename: &str, output_fastq_filename: &str, max_distance: usize) {

        let whitelist = whitelist::read_whitelist(whitelist_filename);
        let barcode_counts = read_barcode_counts(counts_filename);
        let barcode_corrector = BarcodeCorrector::from(&whitelist, &barcode_counts);
        let whitelist: HashSet<&[u8]> = whitelist.iter().map(|s| s.as_bytes()).collect();
    
        let fastq_in = BufReader::new(MultiGzDecoder::new(File::open(input_fastq_filename).unwrap()));
        let fastq_reader = fastq::Reader::from_bufread(fastq_in);
    
        let fastq_out = BufWriter::new(GzEncoder::new(File::create(output_fastq_filename).unwrap(), Compression::fast()));
        let mut fastq_writer = fastq::Writer::from_bufwriter(fastq_out);
    
        let mut matched_whitelist_before_correction: usize = 0;
        let mut matched_whitelist_after_correction: usize = 0;
        let mut total: usize = 0;
    
        for result in fastq_reader.records() {
            total += 1;
    
            let record = result.unwrap();
    
            if whitelist.contains(&record.seq()) {
                matched_whitelist_before_correction += 1;
                matched_whitelist_after_correction += 1;
                let new_description = format!("CR:Z:{}\tCB:Z:{}\tCY:Z:{}", String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(record.qual().to_vec()).unwrap());
    
                fastq_writer.write(record.id(), Some(&new_description), record.seq(), record.qual()).unwrap();
            } else {
                let corrected = barcode_corrector.correct_barcode(record.seq(), record.qual(), max_distance);
    
                let new_description = match corrected {
                    Some(x) => {
                        matched_whitelist_after_correction += 1;
                        format!("CR:Z:{}\tCB:Z:{}\tCY:Z:{}", String::from_utf8(record.seq().to_vec()).unwrap(), x, String::from_utf8(record.qual().to_vec()).unwrap())
                    },
                    None => {
                        format!("CR:Z:{}\tCY:Z:{}", String::from_utf8(record.seq().to_vec()).unwrap(), String::from_utf8(record.qual().to_vec()).unwrap())
                    },
                };
                
                fastq_writer.write(record.id(), Some(&new_description), record.seq(), record.qual()).unwrap();
            }
            
            if total % 1_000_000 == 0 {
                info!("Processed {total} records so far; {matched_whitelist_before_correction} matched whitelist before correction, {matched_whitelist_after_correction} matched whitelist after correction");
            }
        }
    
        fastq_writer.flush().unwrap();
    
    }