#!/usr/bin/env python

import pysam
import random
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('bam_out_prefix', help = 'Output bams will have names {bam_out_prefix}chunk_N.bam')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

def get_barcode(read):
    if not read.has_tag('CB'):
        return 'no_barcode'
    else:
        return read.get_tag('CB')

observed_barcodes = set()
logging.info('Determining all barcodes present.')

with pysam.AlignmentFile(args.bam_in, 'rb') as old_bam:
    count = 0
    for read in old_bam.fetch(until_eof = True):
        count += 1
        observed_barcodes.add(get_barcode(read))
        if count % 1000000 == 0:
            logging.info('Processed {} reads (identified {} barcodes)...'.format(count, len(observed_barcodes)))


# split into chunks of ~10k
all_barcodes = list(observed_barcodes)
random.shuffle(all_barcodes)

chunk = 1
barcode_to_chunk = {}
for count, barcode in enumerate(all_barcodes, 1):
    if count % 10000 == 0:
        chunk += 1
    barcode_to_chunk[barcode] = chunk


logging.info('Creating new bams')
with pysam.AlignmentFile(args.bam_in, 'rb') as old_bam:
    logging.info('Opening {} bam files'.format(len(set(barcode_to_chunk.values()))))
    file_handles = {chunk: pysam.AlignmentFile(f'{args.bam_out_prefix}chunk{chunk}.bam', 'wb', template=old_bam) for chunk in set(barcode_to_chunk.values())}
    count = 0
    for read in old_bam.fetch(until_eof = True):
        count += 1
        if count % 1000000 == 0:
            logging.info('Processed {} reads...'.format(count))
        chunk = barcode_to_chunk[get_barcode(read)]
        file_handles[chunk].write(read)
    for fh in file_handles.values():
        fh.close()
    
logging.info('Done.')