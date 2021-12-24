#!/usr/bin/env python

import os
import sys
import argparse
import logging
import gzip

parser = argparse.ArgumentParser(description='Split a fastq file into N chunks')
parser.add_argument('fastq', help='Fastq file to split. Should be gzipped')
parser.add_argument('n', type=int, help='Number of chunks to create')
parser.add_argument('prefix', help='Prefix for output files (names will be: {PREFIX}chunk_1.fastq, etc')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

with gzip.open(args.fastq, 'rt') as fastq:
    line_count = 0
    record_count = 0
    current_chunk = 1
    
    FHs = dict()
    for chunk in range(1, args.n+1):
        new_file = f'{args.prefix}chunk_{chunk}.fastq'
        logging.info(f'Creating new file: {new_file}')
        FHs[chunk] = open(new_file, 'w')

    for line in fastq:
        line_count += 1
        if (line_count - 1) % 4 == 0:
            record_count += 1
            current_chunk = ((record_count - 1) % args.n) + 1
            if record_count % 1000000 == 0:
                logging.info(f'Processed {record_count} records')
        FHs[current_chunk].write(line)

for fh in FHs.values():
    fh.close()

logging.info('Done')
