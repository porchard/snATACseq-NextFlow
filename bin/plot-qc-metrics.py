#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from skimage.filters import threshold_multiotsu
import numpy as np
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--prefix', default='qc.', help='Prefix for output files (default: "qc.")')
parser.add_argument('ataqv_metrics', help='Path to file of ataqv metrics')
args = parser.parse_args()

ATAQV_METRIC_FILE = args.ataqv_metrics
PREFIX = args.prefix

@ticker.FuncFormatter
def read_count_formatter(x, pos):
    if x >= 1e9:
        return '{}B'.format(x/1e9)
    if x >= 1e6:
        return '{}M'.format(x/1e6)
    if x >= 1e3:
        return '{}k'.format(x/1e3)
    else:
        return x


qc = pd.read_csv(ATAQV_METRIC_FILE, sep='\t', names=['barcode', 'metric', 'value'], header=None).pivot(index='barcode', columns='metric', values='value')

# try to infer HQAA threshold
# do on logscale
values = np.log10(qc[qc.hqaa>100].hqaa).values
values = values.reshape((len(values),1))
thresholds = threshold_multiotsu(image=values, classes=3, nbins=256)
# convert back to linear scale
thresholds = [pow(10, i) for i in thresholds]
HQAA_THRESHOLD = round(thresholds[1])

fig, ax = plt.subplots()
sns.histplot(x='hqaa', data=qc[qc.hqaa>10], ax=ax, log_scale=True)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--', label='HQAA threshold = {:,}'.format(HQAA_THRESHOLD))
ax.legend()
ax.set_xlabel('HQAA')
fig.tight_layout()
fig.savefig(f'{PREFIX}hqaa-distribution.png', dpi=300)
fig.clf()


suggested_thresholds = pd.DataFrame({'metric': ['min_HQAA'], 'threshold': HQAA_THRESHOLD})
suggested_thresholds.to_csv(f'{PREFIX}suggested-thresholds.tsv', sep='\t', index=False)


fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(6*3, 6*2))

# HQAA vs TSS enrichment
ax = axs[0,0]
sns.scatterplot(x='hqaa', y='tss_enrichment', alpha=0.01, data=qc[qc.hqaa>=10], edgecolor=None, ax=ax)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=10)
ax.set_ylim(0.1, 100)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Pass filter reads')
ax.set_ylabel('TSS enrichment')
ax.grid(True)

ax = axs[1,0]
sns.scatterplot(x='hqaa', y='tss_enrichment', alpha=0.01, data=qc[qc.hqaa>=10], edgecolor=None, ax=ax)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(left=1000, right=100000)
ax.set_ylim(0.1, 100)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Pass filter reads')
ax.set_ylabel('TSS enrichment')
ax.grid(True)

# HQAA vs max_fraction_reads_from_single_autosome
ax = axs[0,1]
sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', alpha=0.01, data=qc[qc.hqaa>=10], edgecolor=None, ax=ax)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--')
ax.set_xscale('log')
ax.set_xlim(left=10)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Pass filter reads')
ax.set_ylabel('Max. fraction reads from single autosome')
ax.grid(True)

ax = axs[1,1]
sns.scatterplot(x='hqaa', y='max_fraction_reads_from_single_autosome', alpha=0.01, data=qc[qc.hqaa>=10], edgecolor=None, ax=ax)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--')
ax.set_xscale('log')
ax.set_xlim(left=1000, right=100000)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Pass filter reads')
ax.set_ylabel('Max. fraction reads from single autosome')
ax.grid(True)

# HQAA vs percent mitochondrial
ax = axs[0,2]
sns.scatterplot(x='hqaa', y='percent_mitochondrial', alpha=0.01, data=qc[qc.hqaa>=10], edgecolor=None, ax=ax)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--')
ax.set_xscale('log')
ax.set_xlim(left=10)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Pass filter reads')
ax.set_ylabel('Percent mitochondrial')
ax.grid(True)

ax = axs[1,2]
sns.scatterplot(x='hqaa', y='percent_mitochondrial', alpha=0.01, data=qc[qc.hqaa>=10], edgecolor=None, ax=ax)
ax.axvline(HQAA_THRESHOLD, color='red', ls='--')
ax.set_xscale('log')
ax.set_xlim(left=1000, right=100000)
ax.set_ylim(0, 20)
ax.xaxis.set_major_formatter(read_count_formatter)
ax.set_xlabel('Pass filter reads')
ax.set_ylabel('Percent mitochondrial')
ax.grid(True)

fig.tight_layout()
fig.savefig(f'{PREFIX}metrics.png')
fig.clf()