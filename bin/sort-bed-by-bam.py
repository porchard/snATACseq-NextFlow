#!/usr/bin/env python
# coding: utf-8

import pybedtools as bt
import logging
import argparse
import pysam

parser = argparse.ArgumentParser(description='', add_help = True)
parser.add_argument('bed', type = str,  help = 'BED file to sort')
parser.add_argument('bam', type = str,  help = 'BAM file (should be coordinate sorted).')
parser.add_argument('--use-header', dest = 'use_header', action = 'store_true', default = False, help = 'Infer the bam sort order from the header')
parser.add_argument('--drop-missing', dest = 'drop_missing', action = 'store_true', default = False, help = 'Drop chromosomes missing from the bam file from the bed file.')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

bam = pysam.AlignmentFile(args.bam)
header = bam.header.to_dict() if hasattr(bam, 'to_dict') else bam.header # for python 2 and python 3 compatibility


# determine the sort order...
logging.info('Determining the sort order in the bam file')
sort_order = list()
chromosomes = set()
if args.use_header:
    header = None
    with pysam.AlignmentFile(args.bam) as bam:
        header = bam.header.to_dict() if hasattr(bam, 'to_dict') else bam.header # for python 2 and python 3 compatibility
    sort_order = [i['SN'] for i in header['SQ']]
    chromosomes = set(sort_order)
else:
    with pysam.AlignmentFile(args.bam) as bam:
        for read in bam.fetch(until_eof = True):
            chrom = bam.get_reference_name(read.reference_id)
            if chrom not in chromosomes:
                sort_order.append(chrom)
                chromosomes.add(chrom)

logging.info('Determined the sort order: {}'.format(', '.join(sort_order)))

# sort the bed file
logging.info('Sorting the bed file')
bed = bt.BedTool(args.bed)
for chrom in sort_order:
    chromosome_bed = bed.filter(lambda i: i.chrom == chrom).saveas().sort()
    for i in chromosome_bed:
        print(str(i).rstrip())
        
# are there any chromosomes in the bed file that weren't in the bam file?
missing_chromosomes = bed.filter(lambda i: i.chrom not in chromosomes).saveas()
if len(missing_chromosomes) > 0:
    missing = set()
    if not args.drop_missing:
        for i in missing_chromosomes.sort():
            print(str(i).rstrip())
            missing.add(i.chrom)
        logging.info('Printed chromosomes {} (missing from the bam file) at the end of the bed file.'.format(', '.join(list(missing))))
    else:
        for i in missing_chromosomes.sort():
            missing.add(i.chrom)
        logging.info('Dropped chromosomes {} (missing from the bam file).'.format(', '.join(list(missing))))


logging.info('Done')
