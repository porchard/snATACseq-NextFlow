#!/usr/bin/env python

import os
import sys
import pysam
import logging
import argparse

parser = argparse.ArgumentParser(description='Given a file of barcode corrections and a bam file with uncorrected barcodes, output a new bam file with corrected barcodes.', add_help = True)
parser.add_argument('bam_in', type = str,  help = 'Input bam file. Uncorrected barcodes should be encoded in the read name as: readname_barcode.')
parser.add_argument('corrections', type = str,  help = 'File of barcode corrections (format tsv: readname, old_barcode, new_barcode.')
parser.add_argument('bam_out', type = str,  help = 'Bam file to write.')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

# load in the corrections...
logging.info('Loading barcode corrections from {args.corrections}'.format(**locals()))

corrections = dict()

with open(args.corrections, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        if not len(line) == 3:
            raise ValueError('Expected 3 tab-separated columns in {args.corrections}'.format(**locals()))
        read_name, old_barcode, new_barcode = line
        corrections[read_name] = (old_barcode, new_barcode)

logging.info('Finished loading barcode corrections')
logging.info('Loaded {} barcode corrections'.format(len(corrections)))

logging.info('Correcting barcodes in the input bam file ({args.bam_in})'.format(**locals()))
# write the new file...
with pysam.AlignmentFile(args.bam_in, 'rb') as bam_in:
    with pysam.AlignmentFile(args.bam_out, 'wb', template=bam_in) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            original_read_name, uncorrected_barcode = read.query_name.split('_')
            read.set_tag('CR', uncorrected_barcode)
            if original_read_name in corrections:
                old_barcode, new_barcode = corrections[original_read_name]
                if not old_barcode == uncorrected_barcode:
                    raise ValueError('The old barcode for {} ({}) according to the corrections file does not match the old barcode encoded in the read name ({})'.format(original_read_name, old_barcode, uncorrected_barcode))
                if new_barcode != 'None':
                    read.query_name = '{}_{}'.format(original_read_name, new_barcode)
                    read.set_tag('CB', new_barcode)
                else:
                    continue
            else:
                read.set_tag('CB', uncorrected_barcode)
            bam_out.write(read)

logging.info('Finished writing new bam file ({args.bam_out})'.format(**locals()))
