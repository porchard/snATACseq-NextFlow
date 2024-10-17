#!/usr/bin/env python

import argparse
import sys

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('metrics', help = 'ataqv tabular metrics')
args = parser.parse_args()

def add_metrics(d):
    d['percent_mitochondrial'] = 100 * d.total_mitochondrial_reads.astype(float) / d.total_reads
    d['percent_mitochondrial_duplicate'] = 100 * d.duplicate_mitochondrial_reads.astype(float) / d.total_mitochondrial_reads
    d['percent_autosomal_duplicate'] = np.where(d.total_autosomal_reads > 0, 100 * d.duplicate_autosomal_reads.astype(float) / d.total_autosomal_reads, np.nan)
    d['percent_hqaa'] = 100.0 * d.hqaa.astype(float) / d.total_reads
    d['percent_properly_paired_and_mapped'] = 100.0 * d.properly_paired_and_mapped_reads.astype(float) / d.total_reads
    d['percent_secondary'] = 100.0 * d.secondary_reads.astype(float) / d.total_reads
    d['percent_supplementary'] = 100.0 * d.supplementary_reads.astype(float) / d.total_reads
    d['percent_duplicate'] = 100.0 * d.duplicate_reads.astype(float) / d.total_reads
    d['percent_unmapped'] = 100.0 * d.unmapped_reads.astype(float) / d.total_reads
    d['percent_unmapped_mate'] = 100.0 * d.unmapped_mate_reads.astype(float) / d.total_reads
    d['percent_qcfailed'] = 100.0 * d.qcfailed_reads.astype(float) / d.total_reads
    d['percent_unpaired'] = 100.0 * d.unpaired_reads.astype(float) / d.total_reads
    d['percent_mapq_0'] = 100.0 * d.reads_mapped_with_zero_quality.astype(float) / d.total_reads
    d['percent_rf'] = 100.0 * d.rf_reads.astype(float) / d.total_reads
    d['percent_ff'] = 100.0 * d.ff_reads.astype(float) / d.total_reads
    d['percent_rr'] = 100.0 * d.rr_reads.astype(float) / d.total_reads
    d['percent_autosomal'] = 100.0 * d.total_autosomal_reads.astype(float) / d.total_reads
    d['percent_mate_separate_chromosome'] = 100.0 * d.reads_with_mate_mapped_to_different_reference.astype(float) / d.total_reads
    d['percent_mate_too_distant'] = 100.0 * d.reads_with_mate_too_distant.astype(float) / d.total_reads
    d['percent_improperly_paired'] = 100.0 * d.reads_mapped_and_paired_but_improperly.astype(float) / d.total_reads
    return d


df = pd.read_csv(args.metrics, sep='\t')
df = df[df['name'] != 'no_barcode']
df = add_metrics(df)
df.to_csv(sys.stdout, index=False, sep='\t')