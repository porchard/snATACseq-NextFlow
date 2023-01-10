#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import argparse
import gzip

parser = argparse.ArgumentParser(description='Compare nuclear barcodes to the barcode whitelist to infer whether they need to be reverse complemented', add_help=True)
parser.add_argument('--whitelist', required=True, help='Nuclear barcode whitelist')
parser.add_argument('--fastq', nargs='+', help='Fastq file(s) containing nuclear barcodes')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')



def open_maybe_gzipped(filename):
    """
    Open a possibly gzipped file.
    
    Parameters
    ----------
    filename: str
    The name of the file to open.

    Returns
    -------
    file
    An open file object.
    """
    with open(filename, 'rb') as test_read:
        byte1, byte2 = test_read.read(1), test_read.read(1)
        if byte1 and ord(byte1) == 0x1f and byte2 and ord(byte2) == 0x8b:
            f = gzip.open(filename, mode='rt')
        else:
            f = open(filename, 'rt')
    return f


whitelist = set()
logging.info('Reading whitelist barcodes')

with open(args.whitelist, 'r') as f:
    for line in f:
        line = line.rstrip()
        whitelist.add(line)

logging.info('Finished reading barcodes')

counts = {f: [0, 0] for f in args.fastq} # matching, total


for f in args.fastq:
    logging.info('Processing fastq file {}'.format(f))

    line_count = 0

    with open_maybe_gzipped(f) as fh:
        for line in fh:
            line_count += 1
            if (line_count - 2) % 4 != 0: # skip everything but the sequence itself
                continue
            line = line.rstrip()
            counts[f][1] += 1
            if line in whitelist:
                counts[f][0] += 1


counts_df = pd.DataFrame([counts[f] + [os.path.basename(f)] for f in counts.keys()], columns=['matching', 'total', 'f'])
counts_df['not_matching'] = counts_df.total - counts_df.matching
counts_df['perc_matching'] = 100 * (counts_df.matching / counts_df.total)
counts_df['label'] = ['{}\n({:,}% match whitelist)'.format(f.replace('.transformed-barcode.fastq.gz', ''), round(100*frac, 1)) for f, frac in zip(counts_df.f, counts_df.matching / counts_df.total)]


fig, axs = plt.subplots(ncols=2, figsize=(5*2, 0.5+len(counts_df)))

ax = axs[0]
sns.barplot(x='matching', y='label', data=counts_df, ax=ax)
ax.set_xlabel('# barcodes matching whitelist')
ax.set_ylabel('Readgroup')

ax = axs[1]
sns.barplot(x='perc_matching', y='label', data=counts_df, ax=ax)
ax.set_xlabel('% barcodes matching whitelist')
ax.set_ylabel('Readgroup')
ax.set_xlim(0, 100)

fig.tight_layout()
fig.savefig('barcode-whitelist-matches.png', facecolor='white', dpi=300)
fig.clf()
