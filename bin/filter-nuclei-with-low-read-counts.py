#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import re
import pysam
import logging

parser = argparse.ArgumentParser()
parser.add_argument('--min-reads', dest = 'min_reads', default = 0, type = int, help = 'Minimum number of reads in order to keep nucleus')
parser.add_argument('bam_in', help = 'Bam file to filter')
parser.add_argument('bam_out', help = 'Bam file to write')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


read_counts = dict()
total_reads = 0

with pysam.AlignmentFile(args.bam_in, 'rb') as f:
    for read in f.fetch(until_eof=True):
        total_reads += 1
        if total_reads % 1000000 == 0:
            logging.info('Processed {} reads'.format(total_reads))
        barcode = read.get_tag('CB')
        if barcode not in read_counts:
            read_counts[barcode] = 0
        read_counts[barcode] += 1

keep = set([key for key, value in read_counts.items() if value >= args.min_reads])
reads_kept = sum([value for value in read_counts.values() if value >= args.min_reads])

logging.info('Keeping {} of {} barcodes ({} of {} reads)'.format(len(keep), len(read_counts), reads_kept, total_reads))

with pysam.AlignmentFile(args.bam_in, 'rb') as f:
    with pysam.AlignmentFile(args.bam_out, 'wb', template = f) as f_new:
        for read in f.fetch(until_eof=True):
            barcode = read.get_tag('CB')
            if barcode in keep:
                f_new.write(read)

