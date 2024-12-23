use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Fastq file of (corrected) barcodes
    #[arg(long)]
    barcode_fastq: String,

    /// Read 1 fastq file
    #[arg(long)]
    read_1_fastq: String,

    /// Read 2 fastq file
    #[arg(long)]
    read_2_fastq: String,

    /// Barcode whitelist
    #[arg(long)]
    whitelist_file: String,

    /// Number of chunks to create
    #[arg(long)]
    n_chunks: usize,
}

fn main() {

    let env = env_logger::Env::default().filter_or("MY_LOG_LEVEL", "info");
    env_logger::init_from_env(env);


    let args = Args::parse();
    chunk_by_barcode::run(&args.barcode_fastq, &args.read_1_fastq, &args.read_2_fastq, &args.whitelist_file, args.n_chunks);
}
