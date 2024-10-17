#!/usr/bin/env python
# coding: utf-8


import os
import logging
import argparse
import gzip


import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns


@ticker.FuncFormatter
def read_count_formatter(x, pos):
    """
    Tick label formatting function that converts labels to B/M/k (billions, millions, thousands).

    Usage:
    ax.xaxis.set_major_formatter(read_count_formatter)
    """
    if abs(x) >= 1e9:
        return '{}B'.format(x/1e9)
    elif abs(x) >= 1e6:
        return '{}M'.format(x/1e6)
    elif abs(x) >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


parser = argparse.ArgumentParser(description='Compare nuclear barcodes to the barcode whitelist to infer whether they need to be reverse complemented', add_help=True)
parser.add_argument('fastq', nargs='+', help='Fastq file(s) containing nuclear barcodes. CB:Z and CR:Z tags should be embedded as descriptions (in the read name line). If CB is present, it is assumed to be the corrected barcode and be on the whitelist.')
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


counts = {f: [0, 0, 0] for f in args.fastq} # matching before correction, matching after correction, total


for f in args.fastq:
    logging.info('Processing fastq file {}'.format(f))

    line_count = 0

    with open_maybe_gzipped(f) as fh:
        for line in fh:
            line_count += 1
            if (line_count - 1) % 4 != 0: # skip everything but the read name
                continue
            tags = ' '.join(line.rstrip().split(' ')[1:])
            tags = {i.split(':')[0]: i.split(':')[2] for i in tags.split('\t')}
            if 'CB' in tags:
                if tags['CR'] == tags['CB']:
                    counts[f][0] += 1
                counts[f][1] += 1
            counts[f][2] += 1


counts_df = pd.DataFrame([counts[f] + [os.path.basename(f)] for f in counts.keys()], columns=['matching_before_correction', 'matching_after_correction', 'total', 'f'])

counts_df = counts_df.melt(id_vars=['f', 'total'])
counts_df['perc_matching'] = 100 * (counts_df.value / counts_df.total)
counts_df['label'] = counts_df.f.str.replace('.corrected.fastq.gz', '')

fig, axs = plt.subplots(ncols=2, figsize=(5*2, 0.5+len(counts_df)))

ax = axs[0]
sns.barplot(x='value', y='label', data=counts_df, ax=ax, hue='variable')
ax.set_xlabel('# barcodes matching whitelist')
ax.set_ylabel('Readgroup')
ax.legend().remove()
ax.xaxis.set_major_formatter(read_count_formatter)

ax = axs[1]
sns.barplot(x='perc_matching', y='label', data=counts_df, ax=ax, hue='variable')
ax.set_xlabel('% barcodes matching whitelist')
ax.set_ylabel('Readgroup')
ax.set_xlim(0, 100)
ax.legend(bbox_to_anchor=(1, 1))

fig.tight_layout()
fig.savefig('barcode-whitelist-matches.png', facecolor='white', dpi=300, bbox_inches='tight')
fig.clf()
