#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.chunks = 1

IONICE = 'ionice -c2 -n7'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
            'hg38': (1..22).collect({it -> 'chr' + it}),
            'rn5': (1..20).collect({it -> 'chr' + it}),
            'rn6': (1..20).collect({it -> 'chr' + it}),
            'mm9': (1..19).collect({it -> 'chr' + it}),
            'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
        'hg38': 'human',
        'rn5': 'rat',
        'rn6': 'rat',
        'mm9': 'mouse',
        'mm10': 'mouse']


def make_excluded_regions_arg (genome) {
    return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}


def is_chimeric (library) {
    return get_genome(library).size() > 1
}


def get_bwa_index (genome) {
    params.bwa_index[genome]
}


def get_genome (library) {
    params.libraries[library].genome
}


def get_tss (genome) {
    params.tss[genome]
}


def get_organism (genome) {
    ORGANISMS[genome]
}


def get_chrom_sizes (genome) {
    params.chrom_sizes[genome]
}


def library_to_readgroups (library) {
    params.libraries[library].readgroups.keySet()
}


def library_and_readgroup_to_fastqs (library, readgroup) {
    params.libraries[library].readgroups[readgroup]
}


process fastqc {

    publishDir "${params.results}/fastqc/before-trim", mode: 'rellink', overwrite: true
    time '24h'

    input:
    tuple val(library), val(readgroup), val(read), path(fastq)

    output:
    tuple path(outfile_1), path(outfile_2)

    script:
    outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html')
    outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')

    """
    fastqc $fastq
    """

}

//{args.prefix}chunk_{chunk}.fastq --> ${library}___${readgroup}.${read}.chunk_${chunk}.fastq
process chunk_fastq {

    maxForks 10
    time '24h'

    input:
    tuple val(library), val(readgroup), val(read), path(fastq)

    output:
    tuple val(library), val(readgroup), val(read), path("*.fastq")

    when:
    params.chunks > 1

    """
    chunk-fastq.py $fastq ${params.chunks} ${library}___${readgroup}.${read}.
    """

}

// Trim/reverse complement barcode if necessary. Necessary transformation inferred based on naive comparison of barcode read to barcode whitelist.
process transform_barcode {

    publishDir "${params.results}/transformed-barcodes", mode: 'rellink'
    tag "${library}-${readgroup}"
    memory '10 GB'

    input:
    tuple val(library), val(readgroup), path(fastq)

    output:
    tuple val(library), val(readgroup), path("${library}___${readgroup}.transformed-barcode.fastq.gz")

    """
    ${IONICE} transform-barcode-maybe-gzip.py --check-first 1000000 $fastq ${params['barcode-whitelist']} | gzip -c > ${library}___${readgroup}.transformed-barcode.fastq.gz
    """

}


process make_barcode_corrections {

    publishDir "${params.results}/corrected-barcodes", mode: 'rellink', overwrite: true
    tag "${library}"
    cpus 3
    memory '10 GB'

    input:
    tuple val(library), path(barcode_fastq)

    output:
    tuple val(library), path("${library}.barcode_corrections.txt")

    """
    ${IONICE} correct-barcodes.py --threads 3 ${params['barcode-whitelist']} ${barcode_fastq.join(' ')} > ${library}.barcode_corrections.txt
    """

}


process trim {

    publishDir "${params.results}/trim", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    time '24h'
    tag "${library}-${readgroup}"

    input:
    tuple val(library), val(readgroup), path(fastq_1), path(fastq_2), path(barcode)

    output:
    tuple val(library), val(readgroup), path("${library}-${readgroup}.1.trimmed.fastq.gz"), path("${library}-${readgroup}.2.trimmed.fastq.gz")

    """
    ${IONICE} cta --append-barcode $barcode $fastq_1 $fastq_2 ${library}-${readgroup}.1.trimmed.fastq.gz ${library}-${readgroup}.2.trimmed.fastq.gz
    """

}


process fastqc_post_trim {

    publishDir "${params.results}/fastqc/after-trim", mode: 'rellink', overwrite: true
    time '24h'

    input:
    tuple val(library), val(readgroup), val(read), path(fastq)

    output:
    tuple path(outfile_1), path(outfile_2)

    script:
    outfile_1 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.html')
    outfile_2 = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')

    """
    fastqc $fastq
    """

}


process bwa {

    publishDir "${params.results}/bwa", mode: 'rellink', overwrite: true
    memory '50 GB'
    cpus 12
    errorStrategy 'retry'
    maxRetries 1
    time '48h'
    tag "${library}-${readgroup}-${genome}"

    input:
    tuple val(library), val(genome), val(readgroup), path(fastq_1), path(fastq_2)

    output:
    tuple val(library), val(readgroup), val(genome), path("${library}-${readgroup}-${genome}.bam")

    """
    bwa mem -I 200,200,5000 -M -t 12 ${get_bwa_index(genome)} ${fastq_1} ${fastq_2} | samtools sort -m 1g -@ 11 -O bam -T sort_tmp -o ${library}-${readgroup}-${genome}.bam -
    """

}


process correct_barcodes_in_bam {

    tag "${library}-${readgroup}-${genome}"
    publishDir "${params.results}/bwa-corrected-barcodes", mode: 'rellink', overwrite: true
    memory '75 GB'
    time '24h'

    input:
    tuple val(library), val(readgroup), val(genome), path(bam), path(corrections)

    output:
    tuple val(library), val(genome), path("${library}-${readgroup}-${genome}.corrected.bam")

    """
    correct-barcodes-in-bam.py $bam $corrections ${library}-${readgroup}-${genome}.corrected.bam
    """

}


process merge_readgroups {

    time '24h'
    publishDir "${params.results}/merge", mode: 'rellink', overwrite: true
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path(bams)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.bam")

    """
    samtools merge ${library}-${genome}.bam ${bams.join(' ')}
    """

}


process mark_duplicates {

    publishDir "${params.results}/mark_duplicates", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    time '24h'
    memory '50 GB'
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path("${library}-${genome}.bam")

    output:
    tuple val(library), val(genome), path("${library}-${genome}.md.bam"), path("${library}-${genome}.md.bam.bai")

    """
    java -Xmx40g -Xms40g -jar \$PICARD_JAR MarkDuplicates TMP_DIR=. I=${library}-${genome}.bam O=${library}-${genome}.md.bam READ_ONE_BARCODE_TAG=CB READ_TWO_BARCODE_TAG=CB ASSUME_SORTED=true MAX_RECORDS_IN_RAM=10000000 METRICS_FILE=${library}-${genome}.metrics VALIDATION_STRINGENCY=LENIENT
    samtools index ${library}-${genome}.md.bam
    """

}


process prune {

    publishDir "${params.results}/prune", mode: 'rellink', overwrite: true
    memory '3 GB'
    time '5h'
    errorStrategy 'retry'
    maxRetries 2
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.pruned.bam")

    """
    ${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $md_bam ${AUTOSOMAL_REFERENCES[genome].join(' ')} > ${library}-${genome}.unsorted.bam && samtools sort -m 2G -o ${library}-${genome}.pruned.bam -T bam-sort -O BAM ${library}-${genome}.unsorted.bam
    """

}




process chunk_single_nucleus_bams {

    time '10h'
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.chunk*.bam")

    """
    ${IONICE} chunk-bam-by-barcode.py $md_bam ${library}-${genome}.
    """

}

process ataqv_single_nucleus {

    publishDir "${params.results}/ataqv/single-nucleus/json", mode: 'rellink', overwrite: true
    errorStrategy 'retry'
    maxRetries 1
    memory { 50.GB * task.attempt }
    time '10h'
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), val(chunk), path(md_bam)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.chunk_${chunk}.ataqv.json.gz"), emit: json
    path("${library}-${genome}.chunk_${chunk}.ataqv.out")

    """
    ${IONICE} samtools index $md_bam && ataqv --name ${library}-${genome} --ignore-read-groups --nucleus-barcode-tag CB --metrics-file ${library}-${genome}.chunk_${chunk}.ataqv.json.gz --tss-file ${get_tss(genome)} ${make_excluded_regions_arg(genome)} ${get_organism(genome)} $md_bam > ${library}-${genome}.chunk_${chunk}.ataqv.out
    """

}

process reformat_ataqv {

    memory { 100.GB * task.attempt }
    time '10h'
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path(json)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.txt")

    """
    extractAtaqvMetric.py --files $json > ${library}-${genome}.txt
    """

}


process concat_ataqv {

    publishDir "${params.results}/ataqv/single-nucleus", mode: 'rellink', overwrite: true
    time '10h'
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path("ataqv.*.txt")

    output:
    tuple val(library), val(genome), path("${library}-${genome}.txt")

    """
    cat ataqv.*.txt | cut -f2-4 > ${library}-${genome}.txt
    """

}


workflow {

    libraries = params.libraries.keySet()

    trim_in_inserts = []
    transform_barcode_in = []
    fastqc_in = []
    chunk_fastq_in = []
    no_chunk_fastq_in = []

    for (library in libraries) {
        for (readgroup in library_to_readgroups(library)) {
            fastqs = library_and_readgroup_to_fastqs(library, readgroup)
            first_insert = fastqs['1']
            second_insert = fastqs['2']
            barcode = fastqs['index']
            if (params.chunks > 1) {
                chunk_fastq_in << [library, readgroup, "barcode", file(barcode)]
                chunk_fastq_in << [library, readgroup, "1", file(first_insert)]
                chunk_fastq_in << [library, readgroup, "2", file(second_insert)]
            } else {
                no_chunk_fastq_in << [library, readgroup, "barcode", file(barcode)]
                no_chunk_fastq_in << [library, readgroup, "1", file(first_insert)]
                no_chunk_fastq_in << [library, readgroup, "2", file(second_insert)]
            }
            fastqc_in << [library, readgroup, "1", file(first_insert)]
            fastqc_in << [library, readgroup, "2", file(second_insert)]
            fastqc_in << [library, readgroup, "barcode", file(barcode)]
        }
    }

    fastqc(Channel.from(fastqc_in))

    // handle the chunking
    chunked_out = chunk_fastq(Channel.from(chunk_fastq_in)).transpose().map({it -> [it[0], it[1] + "___" + it[3].getName().tokenize('.')[-2], it[2], it[3]]})
    fastq_in = chunked_out.mix(Channel.from(no_chunk_fastq_in))
    first_insert = fastq_in.filter({it -> it[2] == '1'}).map({it -> [it[0], it[1], it[3]]})
    second_insert = fastq_in.filter({it -> it[2] == '2'}).map({it -> [it[0], it[1], it[3]]})
    barcodes = fastq_in.filter({it -> it[2] == 'barcode'}).map({it -> [it[0], it[1], it[3]]})
    
    transformed_barcodes = transform_barcode(barcodes)
    corrected_barcodes = transformed_barcodes.map({it -> [it[0], it[2]]}).groupTuple(sort: true) | make_barcode_corrections
    trimmed = first_insert.combine(second_insert, by: [0, 1]).combine(transformed_barcodes, by: [0, 1]) | trim
    
    // fastqc on trimmed barcodes
    trimmed.map({it -> [it[0], it[1], ['1', '2'], [it[2], it[3]]]}).transpose() | fastqc_post_trim

    // map
    tmp = []
    for (library in libraries) {
        for(genome in get_genome(library)){
            tmp << [library, genome]
        }
    }

    mapped = Channel.from(tmp).combine(trimmed, by: 0) | bwa

    // processed mapped
    bam_barcodes_corrected = mapped.combine(corrected_barcodes, by: 0) | correct_barcodes_in_bam
    md_bams = bam_barcodes_corrected.groupTuple(by: [0, 1], sort: true) | merge_readgroups | mark_duplicates
    prune(md_bams)
    sn_ataqv = (((md_bams | chunk_single_nucleus_bams).transpose().map({it -> [it[0], it[1], it[2].getName().tokenize('.')[-2].replaceAll('chunk', ''), it[2]]}) | ataqv_single_nucleus).json | reformat_ataqv).groupTuple(by: [0, 1]) | concat_ataqv
}