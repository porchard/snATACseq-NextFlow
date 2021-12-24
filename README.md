# NextFlow pipeline for 10X snATAC-seq data

## Dependencies
If you have Singularity installed, you can use the config provided here ('Singularity') to build a container with all the dependencies.

Otherwise, you'll need to have the following installed:
1. biopython
2. bwa
3. picardtools
4. fastqc
5. samtools (v. 1.10 or above)
6. pysam
7. ataqv (v. 1.3.0 or above)
8. cta (the forked version on the porchard GitHub)

I've used this pipeline with NextFlow v. 20.10.0

## Configuration
Paths to various generic files (e.g., bwa indices) must be included in the nextflow.config file -- check that file and change paths accordingly. These include:

1. Blacklist bed files for each genome
2. Chrom size files for each genome
3. BWA indices
4. TSS files (BED6 files denoting TSS positions)
5. Path to the barcode whitelist (the 10X whitelist is included in this repo)

You'll also need to set the params.results variable -- either in the nextflow.config file itself, or on the command line when you run the pipeline ('--results /path/to/results').

You can split the fastq files into chunks using the --chunks parameter (default: 1, meaning no chunking). In the case of very large fastq files this can speed up processing.

To reduce memory usage of ataqv, we filter out nuclei with low read counts before running ataqv. The minimum read threshold is set in the nextflow.config file.

Lastly, you'll need to include information about each ATAC-seq library, including the genome(s) for the species that each library includes, and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json. Note that for each readgroup, three fastq files are required -- the first and second insert reads ('1' and '2'), and the read with the nuclear barcode ('index')

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -with-singularity /path/to/Singularity.simg -params-file library-config.json --results /path/to/results /path/to/main.nf
```
