use std::fs::File;
use std::collections::HashSet;
use std::io::Read;
use itertools::Itertools;

pub fn read_whitelist(file: &str) -> HashSet<String> {
    let mut whitelist_file = File::open(file).unwrap();
    let mut whitelist: String = String::new();
    whitelist_file.read_to_string(&mut whitelist).unwrap();
    let whitelist: HashSet<String> = whitelist.split("\n").map(|s| s.trim_end().to_string()).collect();
    whitelist
}


pub fn barcode_length(whitelist: &HashSet<String>) -> usize {
    // infer the expected barcode length
    let barcode_lengths: Vec<usize> = whitelist.iter().map(|s| s.len()).dedup().collect();
    assert_eq!(barcode_lengths.len(), 1);
    barcode_lengths[0]
}