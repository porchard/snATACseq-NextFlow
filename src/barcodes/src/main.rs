use clap::{Parser,Subcommand};


#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract barcodes (e.g., 10X cell barcodes) from sequencing reads. 
    /// Depending on the experimental workflow used, barcodes are sometimes embedded in reads 
    /// containing sequence beyond just the barcode, and the barcode may be reverse complemented. 
    /// The `parse-barcodes` option takes a fastq file and a barcode whitelist, and tries to 
    /// infer the location and orientation (reverse complemented or not) of the barcode in the 
    /// read by examining the first 10k records, and write a new fastq file containing the parsed
    /// and properly-oriented barcode.
    ParseBarcodes {
        /// Input fastq file
        #[arg(long)]
        fastq_in: String,

        /// Output fastq file
        #[arg(long)]
        fastq_out: String,

        /// Barcode whitelist
        #[arg(long)]
        whitelist: String,

        /// Output barcode counts file
        #[arg(long)]
        counts: String,
    },
    /// Correct barcodes, using an algorithm similar to that employed in CellRanger's 
    /// ATAC workflow.
    CorrectBarcodes {
        /// Input fastq file
        #[arg(long)]
        fastq_in: String,

        /// Output fastq file
        #[arg(long)]
        fastq_out: String,

        /// Barcode whitelist
        #[arg(long)]
        whitelist: String,

        /// Barcode counts
        #[arg(long)]
        counts: String,

        /// Max Hamming distance
        #[arg(long)]
        max_distance: usize,
    }
}



fn main() {

    let env = env_logger::Env::default().filter_or("MY_LOG_LEVEL", "info");
    env_logger::init_from_env(env);

    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::ParseBarcodes {fastq_in, fastq_out, whitelist, counts}) => {
            barcodes::transform::transform_fastq_file(fastq_in, whitelist, fastq_out, 10000, counts);
        },
        Some(Commands::CorrectBarcodes {fastq_in, fastq_out, whitelist, counts, max_distance}) => {
            barcodes::correct::correct_barcodes_in_fastq(fastq_in, whitelist, counts, fastq_out, *max_distance);
        },
        None => {}
    }
}
