#!/usr/bin/env python
# coding: utf-8


import pybedtools as bt
import gzip
import logging
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('fragment_file')
parser.add_argument('peaks_bed')
parser.add_argument('prefix')
args = parser.parse_args()

FRAGMENT_FILE = args.fragment_file
PEAKS_FILE = args.peaks_bed
PREFIX = args.prefix


barcodes = set()
features = set()
counts = dict() # feature --> barcode --> count

barcodes_chrom_order = []
features_chrom_order = []

# read through the fragment file to create barcodes.tsv
line_count = 0
with gzip.open(FRAGMENT_FILE, 'rt') as ff:
    for line in ff:
        line_count += 1
        if line_count % 1000000 == 0:
            logging.info('Processing line {:,} of fragment file'.format(line_count))
        chrom, start, end, barcode, count = line.rstrip().split()
        if chrom not in barcodes_chrom_order:
            barcodes_chrom_order.append(chrom)
        barcodes.add(barcode)

# read through the peak file to create features.tsv
line_count = 0
with open(PEAKS_FILE, 'rt') as ff:
    for line in ff:
        line_count += 1
        if line_count % 1000000 == 0:
            logging.info('Processing line {:,} of peaks file'.format(line_count))
        chrom, start, end = line.rstrip().split()[:3]
        feature = f'{chrom}_{start}_{end}'
        if chrom not in features_chrom_order:
            features_chrom_order.append(chrom)
        features.add(feature)


# assert(features_chrom_order == barcodes_chrom_order)
barcodes_list = list(sorted(barcodes))
features_list = list(sorted(features))
barcodes = {barcode: i for i, barcode in enumerate(barcodes_list, 1)}
features = {feature: i for i, feature in enumerate(features_list, 1)}

# intersect and write the output
line_count = 0
for i in bt.BedTool(FRAGMENT_FILE).intersect(PEAKS_FILE, sorted=True, wa=True, wb=True, stream=True):
    line_count += 1
    tmp = str(i).rstrip().split()
    barcode_chrom, fragment_start, fragment_end, barcode, _, feature_chrom, feature_start, feature_end = tmp[0:8]
    (fragment_start, fragment_end, feature_start, feature_end) = [int(i) for i in (fragment_start, fragment_end, feature_start, feature_end)]
    feature = f'{feature_chrom}_{feature_start}_{feature_end}'
    # check that one or two cut end overlaps the peak and the peak isn't just in the middle of the fragment
    count = 0
    if fragment_start >= feature_start and fragment_start <= feature_end:
        count += 1
    if fragment_end >= feature_start and fragment_end <= feature_end:
        count += 1
    if count > 0:
        if feature not in counts:
            counts[feature] = dict()
        if barcode not in counts[feature]:
            counts[feature][barcode] = 0
        counts[feature][barcode] += count
        

counts_list = []
for feature, d in counts.items():
    for barcode, count in d.items():
        counts_list.append([str(features[feature]), str(barcodes[barcode]), str(count)])

with open(f'{PREFIX}barcodes.tsv', 'w') as fh:
    for b in barcodes_list:
        fh.write(b + '\n')

with open(f'{PREFIX}features.tsv', 'w') as fh:
    for f in features_list:
        fh.write(f + '\n')

with open(f'{PREFIX}matrix.mtx', 'w') as fh:
    # generate the header
    header = ['%%MatrixMarket matrix coordinate integer general', '%', ' '.join([str(len(features)), str(len(barcodes)), str(len(counts_list))])]
    fh.write('\n'.join(header) + '\n')

    # write the final file
    for count in counts_list:
        fh.write(' '.join(count) + '\n')


# In[ ]:




