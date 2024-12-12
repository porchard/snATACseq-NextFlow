# Barcode preprocessing

This is a small preprocessor that can be used for the following:

1. Extracting barcodes (e.g., 10X cell barcodes) from sequencing reads. Depending on the experimental workflow used, barcodes are sometimes embedded in reads containing sequence beyond just the barcode, and the barcode may be reverse complemented. The `parse-barcodes` option takes a fastq file and a barcode whitelist, and tries to infer the location and orientation (reverse complemented or not) of the barcode in the read, and write a new fastq file containing the parsed and properly-oriented barcode.

2. Correcting barcodes, using the same algorithm used in CellRanger's ATAC workflow (described below; `correct-barcodes` option)

## Installation

You must have Rust installed. You can compile using the following commands:

```
cargo build --release
```

This will create a binary at `./target/release/barcodes`.

## Usage

### Extracting barcodes

The barcode extraction algorithm is as follows:
* Infer the location of the barcode in the read based on the first 10k records. For each record, check each kmer in the record (where k = barcode length, as inferred from the whitelist) against the barcode whitelist. Also check the reverse complement of each kmer against the whitelist. Identify the kmer (and whether or not that kmer must be reverse complemented) that gives the most matches to the whitelist.
* For all of the records in the input fastq file, extract that kmer (and possibly reverse complement) and write it to the output fastq file.

Barcode extraction can be run as follows:

```
./target/release/barcodes parse-barcodes --fastq-in barcodes.fastq.gz --whitelist barcode-whitelist.txt --counts barcode-counts.txt --fastq-out parsed.fastq.gz
```

Where:
* barcodes.fastq.gz is the fastq file or barcode reads. Must be gzipped.
* barcode-whitelist.txt is the barcode whitelist (e.g., a 10X barcode whitelist). Must be unzipped.
* barcode-counts.txt is the path to a two-column TSV file (no header) that will be created, listing barcode in the first column and the number of times that barcode is observed (after extracting barcodes from the reads) in the second column. Must be unzipped. Includes all barcodes that are observed at least once, regardless of whether or not it is on the whitelist.
* parsed.fastq.gz is the output fastq file, containing only the barcode sequences.

### Correcting barcodes

The barcode correction algorithm is as follows, for each barcode not in the whitelist:
* Find all whitelisted barcodes within Hamming distance N (N is user-defined). These are the potential corrections.
* For each potential correction:
  * Use the uncorrected barcode phred scores to calculate the probability of errors at the positions differing between the potential correction and the uncorrected barcode (taking the product of these probabilities if there are > 1 position differing)
  * Multiply the probability calculated above by the number of times that potential correction appears in the dataset (i.e. this is similar to a prior; the more a whitelisted barcode appears in the data, the higher the probability that a given uncorrected barcode will be derived from this whitelisted barcode. A pseudocount of 1 is added to each of the potential corrections).
* These values are then normalized such that they sum to 1 across the potential corrections for a non-whitelisted barcode (i.e., x / sum(x)). If any of the values are >= 0.975, the non-whitelisted barcode is corrected to the corresponding potential correction.

Barcode correction can be run as follows:

```
./target/release/barcodes correct-barcodes --fastq-in barcodes.fastq.gz --whitelist barcode-whitelist.txt --counts barcode-counts.txt --fastq-out corrected.fastq.gz --max-distance 2
```

Where:
* barcodes.fastq.gz is the fastq file or barcode reads (each read should be a properly oriented barcode, without any extra sequence). Must be gzipped.
* barcode-whitelist.txt is the barcode whitelist (e.g., a 10X barcode whitelist). Must be unzipped.
* barcode-counts.txt is a two-column TSV file (no header) listing barcode in the first column and the number of times that barcode is observed in the second column. Must be unzipped.
* --max-distance is the maximum Hamming distance between an uncorrected barcode and the whitelisted barcodes that should be considered for the correction.
* corrected.fastq.gz is the fastq file including barcode corrections.