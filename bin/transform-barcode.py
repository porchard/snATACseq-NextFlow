#!/usr/bin/env python
# coding: utf-8

import os
import sys
import gzip
import logging
import argparse

parser = argparse.ArgumentParser(description='Compare nuclear barcodes to the barcode whitelist to infer whether they need to be reverse complemented', add_help=True)
parser.add_argument('--check-first', default=1000000, type=int, help='Check no more than this number of barcodes (default: 1000000)')
parser.add_argument('fastq', help='Fastq file containing nuclear barcodes (gzipped)')
parser.add_argument('whitelist', help='Nuclear barcode whitelist')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

COMPLEMENTS = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N'
}


def reverse_complement(s):
    return ''.join([COMPLEMENTS[x.upper()] for x in s][::-1])


whitelist_barcodes = set()
whitelist_barcodes_rc = set()
logging.info('Reading whitelist barcodes and determining their reverse complements')

with open(args.whitelist, 'r') as f:
    for line in f:
        line = line.rstrip()
        whitelist_barcodes.add(line)
        whitelist_barcodes_rc.add(reverse_complement(line))

logging.info('Finished reading barcodes')

# determine barcode length
barcode_lengths = set()
for i in whitelist_barcodes:
    barcode_lengths.add(len(i))
if not len(barcode_lengths) == 1:
    raise ValueError('Barcodes are not all of uniform length')
barcode_length = list(barcode_lengths)[0]

logging.info('Processing fastq file {}'.format(args.fastq))

match_counts = dict() # offset --> [number_in_whitelist, number_in_whitelist_rc]
line_count = 0
record_count = 0
barcode_counts = dict()
with gzip.open(args.fastq, 'rt') as f:
    for line in f:
        line_count += 1
        if (line_count - 2) % 4 != 0: # skip everything but the sequence itself
            continue
        record_count += 1
        line = line.rstrip()
        if line not in barcode_counts:
            barcode_counts[line] = 0
        barcode_counts[line] += 1
        if len(line) < barcode_length:
            raise ValueError('Barcode reads are shorted than the barcodes in the whitelist')
        for offset in range(0, len(line) - barcode_length + 1):
            if offset not in match_counts:
                match_counts[offset] = [0, 0]
            potential_barcode = line[offset:(offset+barcode_length)]
            if potential_barcode in whitelist_barcodes:
                match_counts[offset][0] += 1
            if potential_barcode in whitelist_barcodes_rc:
                match_counts[offset][1] += 1
        if record_count == args.check_first:
            break
        if record_count % 1000000 == 0:
            for offset in match_counts:
                in_whitelist_count = match_counts[offset][0]
                in_whitelist_rc_count = match_counts[offset][1]
                percent_in_whitelist = round(100*float(in_whitelist_count) / record_count, 2)
                percent_in_whitelist_rc = round(100*float(in_whitelist_rc_count) / record_count, 2)
                logging.info('Processed {} records so far. For offset {}, {} ({}%) in whitelist, {} ({}%) in whitelist reverse complemented'.format(record_count, offset, in_whitelist_count,  percent_in_whitelist, in_whitelist_rc_count, percent_in_whitelist_rc))

logging.info('Finished examining records')
logging.info('Determining correct transformation.')
OFFSET = 0
RC = False
MAX_MATCH_COUNT = 0
for offset in match_counts:
    matches, rc_matches = match_counts[offset]
    if matches > MAX_MATCH_COUNT:
        MAX_MATCH_COUNT = matches
        OFFSET = offset
        RC = False
    if rc_matches > MAX_MATCH_COUNT:
        MAX_MATCH_COUNT = rc_matches
        OFFSET = offset
        RC = True
logging.info('Determined that barcode begins at position {} and {} reverse complemented'.format(OFFSET, 'is' if RC else 'is not'))

logging.info('Transforming barcode')
line_count = 0
record_count = 0
match_whitelist_count = 0
with gzip.open(args.fastq, 'rt') as f:
    for line in f:
        line = line.rstrip()
        line_count += 1
        if line_count % 4 == 1:
            # name
            print(line)
        if line_count % 4 == 2:
            # read
            if len(line) > barcode_length:
                # trim
                line = line[OFFSET:(OFFSET+barcode_length)]
            if RC:
                line = reverse_complement(line)
            record_count += 1
            match_whitelist_count += 1 if line in whitelist_barcodes else 0
            print(line)
        if line_count % 4 == 3:
            # '+'
            print(line)
        if line_count % 4 == 0:
            # phred
            if len(line) > barcode_length:
                # trim
                line = line[OFFSET:(OFFSET+barcode_length)]
            if RC:
                line = line[::-1]
            print(line)

logging.info('Finished transforming barcodes')
logging.info('{}% of transformed barcodes matched the whitelist'.format(round(100*match_whitelist_count / record_count, 2)))
