#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.chunks = 1

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
    'hg38': (1..22).collect({it -> 'chr' + it}),
    'rn4': (1..20).collect({it -> 'chr' + it}),
    'rn5': (1..20).collect({it -> 'chr' + it}),
    'rn6': (1..20).collect({it -> 'chr' + it}),
    'rn7': (1..20).collect({it -> 'chr' + it}),
    'mm9': (1..19).collect({it -> 'chr' + it}),
    'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
    'hg38': 'human',
    'rn4': 'rat',
    'rn5': 'rat',
    'rn6': 'rat',
    'rn7': 'rat',
    'mm9': 'mouse',
    'mm10': 'mouse'
]

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'rn7': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]


def has_blacklist (genome) {
    return params.blacklist.containsKey(genome)
}


def get_blacklists (genome) {
    if (params.blacklist[genome] instanceof String) {
        return [params.blacklist[genome]]
    } else {
        return params.blacklist[genome]
    }
}


def make_excluded_regions_arg (genome) {
    return get_blacklists(genome).collect({'--excluded-region-file ' + it}).join(' ')
}


def is_chimeric (library) {
    return get_genome(library).size() > 1
}


def get_bwa_index (genome) {
    return params.bwa_index[genome]
}


def get_genome (library) {
    return (params.libraries[library].genome instanceof List) ? params.libraries[library].genome : [params.libraries[library].genome]
}


def get_tss (genome) {
    return params.tss[genome]
}


def get_organism (genome) {
    return ORGANISMS[genome]
}


def get_chrom_sizes (genome) {
    return params.chrom_sizes[genome]
}


def get_macs2_genome_size (genome) {
    return MACS2_GENOME_SIZE[genome]
}


def library_to_readgroups (library) {
    return params.libraries[library].readgroups.keySet()
}


def library_and_readgroup_to_fastqs (library, readgroup) {
    return params.libraries[library].readgroups[readgroup]
}


process fastqc {

    publishDir "${params.results}/fastqc/before-trim"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
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


process multiqc {

    publishDir "${params.results}/multiqc/before-trim"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
    time '5h'

    input:
    path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


// Trim/reverse complement barcode if necessary. Necessary transformation inferred based on naive comparison of barcode read to barcode whitelist.
process transform_barcode {

    tag "${library}-${readgroup}"
    cpus 1
    time '24h'
    memory '5 GB'
    container 'docker://porchard/snatacseq-nextflow-barcodes:20250702'

    input:
    tuple val(library), val(readgroup), path(fastq), path(whitelist)

    output:
    tuple val(library), val(readgroup), path("${library}___${readgroup}.transformed-barcode.fastq.gz"), emit: transformed
    tuple val(library), path("${library}___${readgroup}.counts.txt"), emit: counts

    """
    barcodes parse-barcodes --fastq-in $fastq --whitelist $whitelist --counts ${library}___${readgroup}.counts.txt --fastq-out ${library}___${readgroup}.transformed-barcode.fastq.gz
    """

}


process correct_barcodes {

    tag "${library}"
    memory '10 GB'
    time '24h'
    container 'docker://porchard/snatacseq-nextflow-barcodes:20250702'
    label 'largemem'

    input:
    tuple val(library), val(readgroup), path(barcode_fastq), path("counts/?_counts.txt"), path(whitelist)

    output:
    tuple val(library), val(readgroup), path("${library}___${readgroup}.corrected.fastq.gz")

    """
    cat counts/* > counts.txt
    barcodes correct-barcodes --fastq-in $barcode_fastq --whitelist $whitelist --counts counts.txt --fastq-out ${library}___${readgroup}.corrected.fastq.gz --max-distance 1
    """

}


process plot_whitelist_matching {

    publishDir "${params.results}/plot-barcodes-matching-whitelist"
    memory '10 GB'
    cpus 1
    time '24h'
    container 'library://porchard/default/general:20220107'
    errorStrategy 'ignore'

    input:
    path(fastq)

    output:
    path("barcode-whitelist-matches.png")

    """
    plot-barcodes-matching-whitelist.py ${fastq.join(' ')}
    """

}


process chunk_fastq {

    container 'docker://porchard/snatacseq-nextflow-chunk-by-barcode:20241223'
    memory '4 GB'
    cpus 1
    time '24h'

    input:
    tuple val(library), val(readgroup), path(read_1), path(read_2), path(barcodes), path(whitelist)

    output:
    tuple val(library), val(readgroup), path("*.fastq.gz")

    """
    chunk-by-barcode --n-chunks ${params.chunks} --read-1-fastq ${read_1} --read-2-fastq ${read_2} --barcode-fastq ${barcodes} --whitelist-file ${whitelist}
    """

}


process trim_unchunked {

    errorStrategy 'retry'
    maxRetries 1
    tag "${library}-${readgroup}-${chunk}"
    container 'docker://porchard/cta:c97ac86'
    memory '4 GB'
    cpus 1
    time '24h'

    input:
    tuple val(library), val(readgroup), path(fastq_1), path(fastq_2), path(barcode)

    output:
    tuple val(library), val(readgroup), val("chunk_1"), path("${library}-${readgroup}-chunk_1.1.trimmed.fastq.gz"), path("${library}-${readgroup}-chunk_1.2.trimmed.fastq.gz")

    """
    cta --strip-description --copy-description $barcode $fastq_1 $fastq_2 ${library}-${readgroup}-chunk_1.1.trimmed.fastq.gz ${library}-${readgroup}-chunk_1.2.trimmed.fastq.gz
    """

}


process trim_chunked {

    errorStrategy 'retry'
    maxRetries 1
    tag "${library}-${readgroup}-${chunk}"
    container 'docker://porchard/cta:c97ac86'
    memory '4 GB'
    cpus 1
    time '24h'

    input:
    tuple val(library), val(readgroup), val(chunk), path(fastq_1), path(fastq_2)

    output:
    tuple val(library), val(readgroup), val(chunk), path("${library}-${readgroup}-${chunk}.1.trimmed.fastq.gz"), path("${library}-${readgroup}-${chunk}.2.trimmed.fastq.gz")

    """
    cta $fastq_1 $fastq_2 ${library}-${readgroup}-${chunk}.1.trimmed.fastq.gz ${library}-${readgroup}-${chunk}.2.trimmed.fastq.gz
    """

}


process fastqc_post_trim {

    publishDir "${params.results}/fastqc/after-trim"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
    time '24h'

    input:
    tuple val(library), val(readgroup), val(read), path(fastq)

    output:
    tuple path("${library}-${readgroup}.${read}_fastqc.html"), path("${library}-${readgroup}.${read}_fastqc.zip")

    """
    zcat ${fastq.join(' ')} | gzip -c > ${library}-${readgroup}.${read}.fastq.gz
    fastqc ${library}-${readgroup}.${read}.fastq.gz
    rm ${library}-${readgroup}.${read}.fastq.gz
    """

}


process multiqc_post_trim {

    publishDir "${params.results}/multiqc/after-trim"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1
    time '5h'

    input:
    path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


process bwa {

    memory '50 GB'
    cpus 12
    errorStrategy 'retry'
    maxRetries 1
    time '72h'
    tag "${library}-${readgroup}-${chunk}-${genome}"
    container 'library://porchard/default/bwa:0.7.15'
    label 'largemem'

    input:
    tuple val(library), val(genome), val(readgroup), val(chunk), path(fastq_1), path(fastq_2)

    output:
    tuple val(library), val(genome), val(chunk), path("${library}-${readgroup}-${chunk}-${genome}.bam")

    """
    bwa mem -C -I 200,200,5000 -M -t 12 ${get_bwa_index(genome)} ${fastq_1} ${fastq_2} | samtools sort -m 1g -@ 11 -O bam -T sort_tmp -o ${library}-${readgroup}-${chunk}-${genome}.bam -
    """

}


process merge_readgroups {

    time '24h'
    tag "${library} ${genome} ${chunk}"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1

    input:
    tuple val(library), val(genome), val(chunk), path(bams)

    output:
    tuple val(library), val(genome), val(chunk), path("${library}-${genome}-${chunk}.bam")

    script:
    cmd = (bams instanceof List && bams.size() > 1) ? "samtools merge ${library}-${genome}-${chunk}.bam ${bams.join(' ')}" : "ln -s ${bams[0]} ${library}-${genome}-${chunk}.bam"

    """
    $cmd
    """

}


process mark_duplicates {

    errorStrategy 'retry'
    maxRetries 1
    time '24h'
    memory '50 GB'
    tag "${library} ${genome} ${chunk}"
    container 'library://porchard/default/general:20220107'
    cpus 1
    label 'largemem'

    input:
    tuple val(library), val(genome), val(chunk), path(bam)

    output:
    tuple val(library), val(genome), val(chunk), path("${library}-${genome}-${chunk}.md.bam"), path("${library}-${genome}-${chunk}.md.bam.bai")

    """
    java -Xmx40g -Xms40g -jar \$PICARD_JAR MarkDuplicates TMP_DIR=. I=$bam O=${library}-${genome}-${chunk}.md.bam READ_ONE_BARCODE_TAG=CB READ_TWO_BARCODE_TAG=CB ASSUME_SORTED=true MAX_RECORDS_IN_RAM=10000000 METRICS_FILE=${library}-${genome}.metrics VALIDATION_STRINGENCY=LENIENT
    samtools index ${library}-${genome}-${chunk}.md.bam
    """

}


process merge_chunked_marked_duplicates {

    publishDir "${params.results}/mark_duplicates"
    time '24h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1

    input:
    tuple val(library), val(genome), val(chunks), path(bams), path(indices)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.bam"), path("${library}-${genome}.bam.bai")

    script:
    cmd = (bams instanceof List && bams.size() > 1) ? "samtools merge ${library}-${genome}.bam ${bams.join(' ')}; samtools index ${library}-${genome}.bam" : "ln -s ${bams[0]} ${library}-${genome}.bam; ln -s ${indices[0]} ${library}-${genome}.bam.bai"

    """
    $cmd
    """

}


process prune {

    publishDir "${params.results}/prune"
    memory '3 GB'
    time '24h'
    errorStrategy 'retry'
    maxRetries 2
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    cpus 1

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.pruned.bam")

    """
    samtools view -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $md_bam ${AUTOSOMAL_REFERENCES[genome].join(' ')} | samtools sort -m 2G -o ${library}-${genome}.pruned.bam -T bam-sort -O BAM
    """

}

process index_pruned {

    publishDir "${params.results}/prune"
    memory '3 GB'
    time '24h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    cpus 1

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple val(library), val(genome), path(bam), path("${bam}.bai")

    """
    samtools index $bam
    """

}


process make_fragment_file {

    publishDir "${params.results}/fragment-file"
    memory '60 GB'
    time '24h'
    tag "${library} ${genome}"
    container "docker://porchard/sinto:20230623"
    cpus 10

    input:
    tuple val(library), val(genome), path(bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.frags.bed.gz"), path("${library}-${genome}.frags.bed.gz.tbi")

    """
    sinto fragments -b $bam -f unsorted.frags.bed --collapse_within --nproc 10 --chunksize 5000000
    sort --parallel=10 -S 25G -k 1,1 -k2,2n unsorted.frags.bed > ${library}-${genome}.frags.bed
    bgzip ${library}-${genome}.frags.bed
    tabix -p bed ${library}-${genome}.frags.bed.gz
    rm unsorted.frags.bed
    """

}

process bamtobed {

    time '24h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1

    input:
    tuple val(library), val(genome), path(bam)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.bed")

    """
    bedtools bamtobed -i $bam > ${library}-${genome}.bed
    """

}


process macs2 {

    publishDir "${params.results}/macs2"
    time '24h'
    tag "${library} ${genome}"
    memory { 25.GB * task.attempt }
    maxRetries 2
    errorStrategy 'retry'
    container 'library://porchard/default/general:20220107'
    cpus 1
    label 'largemem'

    input:
    tuple val(library), val(genome), path(bed)

    output:
    tuple val(library), val(genome), path("${library}-${genome}_peaks.broadPeak"), emit: peaks
    tuple val(library), val(genome), path("${library}-${genome}_treat_pileup.bdg"), emit: bedgraph

    """
    macs2 callpeak -t $bed --outdir . --SPMR -f BED -n ${library}-${genome} -g ${get_macs2_genome_size(genome)} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
    """

}

process blacklist_filter_peaks {

    publishDir "${params.results}/macs2"
    time '3h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    cpus 1

    input:
    tuple val(library), val(genome), path(peaks)

    output:
    tuple val(library), val(genome), path("${library}-${genome}_peaks.broadPeak.noblacklist")

    when:
    has_blacklist(genome)

    """
    bedtools intersect -a $peaks -b ${get_blacklists(genome).join(' ')} -v > ${library}-${genome}_peaks.broadPeak.noblacklist
    """

}


process bigwig {

    time '24h'
    publishDir "${params.results}/bigwig"
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    memory { 20.GB * task.attempt }
    maxRetries 2
    errorStrategy 'retry'
    cpus 1
    label 'largemem'

    input:
    tuple val(library), val(genome), path(bedgraph)

    output:
    tuple val(genome), path("${library}-${genome}.bw")

    """
    LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph > sorted.bedgraph
    bedClip sorted.bedgraph ${get_chrom_sizes(genome)} clipped.bedgraph
    bedGraphToBigWig clipped.bedgraph ${get_chrom_sizes(genome)} ${library}-${genome}.bw
    rm sorted.bedgraph clipped.bedgraph
    """

}


process plot_signal_at_tss {

    publishDir "${params.results}/bigwig/plot"
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    tag "${genome}"
    container 'library://porchard/default/general:20220107'
    cpus 1
    time '24h'

    input:
    tuple val(genome), path(bw)

    output:
    path("*.png") optional true

    """
    plot-signal-at-tss.py --genes ${params.plot_signal_at_genes.join(' ')} --tss-file ${get_tss(genome)} --bigwigs ${bw.join(' ')}
    """

}

process ataqv_single_nucleus {

    errorStrategy 'retry'
    maxRetries 1
    memory { 5.GB * task.attempt }
    time '24h'
    tag "${library} ${genome} ${chunk}"
    container 'docker://porchard/ataqv:1.5.0'

    input:
    tuple val(library), val(genome), val(chunk), path(md_bam), path(bam_index)

    output:
    tuple val(library), val(genome), path("${library}-${genome}-${chunk}.ataqv.txt.gz")

    """
    ataqv --name ${library}-${genome} --ignore-read-groups --nucleus-barcode-tag CB --metrics-file ${library}-${genome}-${chunk}.ataqv.txt.gz --tss-file ${get_tss(genome)} ${make_excluded_regions_arg(genome)} ${get_organism(genome)} $md_bam > ${library}-${genome}-${chunk}.ataqv.out
    """

}


process merge_chunked_ataqv_single_nucleus {

    publishDir "${params.results}/ataqv/single-nucleus"
    errorStrategy 'retry'
    maxRetries 1
    time '1h'
    tag "${library} ${genome}"

    input:
    tuple val(library), val(genome), path(metrics)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.ataqv.txt.gz")

    """
    zcat ${metrics[0]} | awk 'NR==1' > header.txt
    zcat ${metrics.join(' ')} | grep -v -f header.txt | cat header.txt - | gzip > ${library}-${genome}.ataqv.txt.gz
    """

}


process add_qc_metrics {

    publishDir "${params.results}/ataqv/single-nucleus"
    time '1h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    memory "7 GB"

    input:
    tuple val(library), val(genome), path(metrics)

    output:
    tuple val(library), val(genome), path("${library}-${genome}.txt")

    """
    add-metrics.py $metrics > ${library}-${genome}.txt
    """

}


process plot_qc_metrics {

    publishDir "${params.results}/ataqv/single-nucleus"
    time '10h'
    tag "${library} ${genome}"
    container 'library://porchard/default/dropkick:20220225'
    memory { 10.GB * task.attempt }
    maxRetries 1
    errorStrategy 'retry'
    cpus 1

    input:
    tuple val(library), val(genome), path(metrics)

    output:
    tuple val(library), val(genome), path("*.png")
    path("*.tsv")

    """
    plot-qc-metrics.py --prefix ${library}-${genome}. $metrics
    """

}


process ataqv_bulk {

    publishDir "${params.results}/ataqv/bulk"
    errorStrategy 'retry'
    maxRetries 1
    memory { 5.GB * task.attempt }
    time '24h'
    tag "${library} ${genome}"
    container 'docker://porchard/ataqv:1.5.0'
    cpus 1

    input:
    tuple val(library), val(genome), path(md_bam), path(bam_index), path(peaks)

    output:
    tuple val(genome), path("${library}-${genome}.ataqv.json.gz"), emit: json
    path("${library}-${genome}.ataqv.out")

    """
    ataqv --name ${library}-${genome} --peak-file $peaks --ignore-read-groups --metrics-file ${library}-${genome}.ataqv.json.gz --tss-file ${get_tss(genome)} ${make_excluded_regions_arg(genome)} ${get_organism(genome)} $md_bam > ${library}-${genome}.ataqv.out
    """

}


process ataqv_bulk_viewer {

    publishDir "${params.results}/ataqv/bulk"
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    time '4h'
    tag "${genome}"
    container 'docker://porchard/ataqv:1.4.0'
    cpus 1

    input:
    tuple val(genome), path(json)

    output:
    path("ataqv-viewer-${genome}")

    """
    mkarv ataqv-viewer-${genome} ${json.join(' ')}
    """

}


process get_peak_counts {

    publishDir "${params.results}/counts"
    memory { 75.GB * task.attempt }
    time '48h'
    tag "${library} ${genome}"
    container 'library://porchard/default/general:20220107'
    cpus 1
    label 'largemem'
    maxRetries 2

    input:
    tuple val(library), val(genome), path(fragments), path(fragments_index), path(peaks)

    output:
    path("${library}-${genome}.atac.*")

    """
    sort -k 1,1 -k2,2n $peaks > peaks.sorted.bed
    fragment-file-to-peak-counts.py $fragments $peaks ${library}-${genome}.atac.
    """

}


workflow {

    libraries = params.libraries.keySet()

    first_insert = []
    second_insert = []
    barcodes = []

    for (library in libraries) {
        for (readgroup in library_to_readgroups(library)) {
            fastqs = library_and_readgroup_to_fastqs(library, readgroup)
            
            first_insert << [library, readgroup, file(fastqs['1'])]
            second_insert << [library, readgroup, file(fastqs['2'])]
            barcodes << [library, readgroup, file(fastqs['index'])]
        }
    }

    first_insert = Channel.from(first_insert)
    second_insert = Channel.from(second_insert)
    barcodes = Channel.from(barcodes)
    whitelist = Channel.fromPath(params['barcode-whitelist'])

    fastqc(first_insert.map({it -> [it[0], it[1], "1", it[2]]}).mix(second_insert.map({it -> [it[0], it[1], "2", it[2]]})).mix(barcodes.map({it -> [it[0], it[1], "barcode", it[2]]}))).flatten().toSortedList() | multiqc

    transformed_barcodes = transform_barcode(barcodes.combine(whitelist))
    corrected_barcodes = transformed_barcodes.transformed.combine(transformed_barcodes.counts.groupTuple(by: 0), by: 0).combine(whitelist) | correct_barcodes
    plot_whitelist_matching(corrected_barcodes.map({it -> it[2]}).toSortedList())


    if (params.chunks > 1) {
        chunked = first_insert.combine(second_insert, by: [0, 1]).combine(corrected_barcodes, by: [0, 1]).combine(whitelist) | chunk_fastq
        chunked = chunked.transpose()
        chunked = chunked.map({it -> [it[0], it[1], it[2].getName().tokenize('.')[0], it[2].getName().tokenize('.')[1], it[2]]})
        chunked_read_1 = chunked.filter({it -> it[2] == 'R1'}).map({it -> [it[0], it[1], it[3], it[4]]})
        chunked_read_2 = chunked.filter({it -> it[2] == 'R2'}).map({it -> [it[0], it[1], it[3], it[4]]})
        trimmed = chunked_read_1.combine(chunked_read_2, by: [0, 1, 2]) | trim_chunked
    } else {
        trimmed = first_insert.combine(second_insert, by: [0, 1]).combine(corrected_barcodes, by: [0, 1]) | trim_unchunked
    }

    
    (trimmed.map({it -> [it[0], it[1], ['1', '2'], [it[3], it[3]]]}).transpose().groupTuple(by: [0, 1, 2]) | fastqc_post_trim).flatten().toSortedList() | multiqc_post_trim

    // map
    tmp = []
    for (library in libraries) {
        for(genome in get_genome(library)){
            tmp << [library, genome]
        }
    }

    mapped = Channel.from(tmp).combine(trimmed, by: 0) | bwa

    // processed mapped
    chunked_md_bams = mapped.groupTuple(by: [0, 1, 2], sort: true) | merge_readgroups | mark_duplicates
    md_bams = merge_chunked_marked_duplicates(chunked_md_bams.groupTuple(by: [0, 1]))
    pruned = prune(md_bams)
    fragment_file = index_pruned(pruned) | make_fragment_file
    peak_calling = bamtobed(pruned) | macs2
    blacklist_filtered_peaks = blacklist_filter_peaks(peak_calling.peaks)
    bigwig(peak_calling.bedgraph).groupTuple() | plot_signal_at_tss
    fragment_file.combine(blacklist_filtered_peaks, by: [0, 1]) | get_peak_counts

    sn_ataqv = (chunked_md_bams | ataqv_single_nucleus).groupTuple(by: [0, 1]) | merge_chunked_ataqv_single_nucleus | add_qc_metrics | plot_qc_metrics
    ataqv_bulk(md_bams.combine(peak_calling.peaks, by: [0, 1])).json.groupTuple() | ataqv_bulk_viewer
}
