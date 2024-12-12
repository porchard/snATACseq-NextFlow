# NextFlow pipeline for 10X snATAC-seq data

## Dependencies
Singularity (v. 3) and NextFlow (>= v. 20.10.0). Containers with the software for each step are pulled from the Sylabs cloud library (https://cloud.sylabs.io/library).


## Configuration
Paths to various generic files (e.g., bwa indices) must be included in the nextflow.config file -- check that file and change paths accordingly. These include:

1. Blacklist bed files for each genome
2. Chrom size files for each genome
3. bwa indices (compatible with bwa v. 0.7.15)
4. TSS files (BED6 files denoting TSS positions)

When launching the pipeline, as shown in the `nextflow` command below, you'll also need to set the following:

1. The location of the results directory (e.g., `--results /path/to/results`)
2. The location of the barcode whitelist (e.g., `--barcode-whitelist /path/to/737K-arc-v1.txt`). The 10X ATAC v1 (737K-cratac-v1.txt) and 10X ATAC multiome (737K-arc-v1.txt) whitelists are included in this repo.

You can split the fastq files into chunks using the --chunks parameter (default: 1, meaning no chunking). In the case of very large fastq files this can speed up processing.

You can generate output plots of the (pseudobulk ATAC) signal at gene TSS by adding gene names to the params.plot_signal_at_genes variable (these gene names must be present in the TSS files). By default only the signal at the GAPDH TSS is plotted.

Lastly, you'll need to include information about each ATAC-seq library, including the genome(s) for the species that each library includes, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json. Note that for each readgroup, three fastq files are required -- the first and second insert reads ('1' and '2'), and the read with the nuclear barcode ('index')

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -resume -params-file library-config.json --barcode-whitelist /path/to/737K-arc-v1.txt --results /path/to/results /path/to/main.nf
```

## Output
* `ataqv/bulk/*.{json.gz,out}`: Pseudobulk ataqv output for each library
* `ataqv/bulk/ataqv-viewer-{genome}`: Pseudobulk ataqv HTML reports
* `ataqv/single-nucleus/*.png`: Plots of per-barcode ataqv metrics
* `ataqv/single-nucleus/*.txt`: Per barcode ataqv metrics in txt format
* `ataqv/single-nucleus/*.suggested-thresholds.tsv`: Suggested min HQAA threshold for the library, based on multi-otsu thresholding of the HQAA distribution
* `bigwig/*.bw`: Pseudobulk bigwig files for each library
* `bigwig/plot/*.png`: Pseudobulk ATAC signal at gene TSS for selected genes
* `macs2/*`: Pseudobulk peak calling results for each library
* `multiqc/*`: multiqc summaries of fastqc results, before and after adapter trimming
* `prune/*`: Filtered bam files
* `plot-barcodes-matching-whitelist`: Plot displaying percentage of barcodes matching barcode whitelist