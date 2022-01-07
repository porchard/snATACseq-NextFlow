#!/usr/bin/env python
# coding: utf-8

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import pygenometracks.tracks as pygtk
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('--genes', required=True, nargs='+', help='List of gene names.')
parser.add_argument('--tss-file', dest='tss_file', required=True, help='TSS bed file.')
parser.add_argument('--bigwigs', required=True, nargs='+', help='Bigwig files')
args = parser.parse_args()

GENES = args.genes
TSS_FILE = args.tss_file
BIGWIGS = args.bigwigs
FLANK = 10e3

logging.info('Loading TSS')
tss = pd.read_csv(TSS_FILE, sep='\t', header=None).iloc[:,0:4]
tss.columns = ['chrom', 'start', 'end', 'gene']

for gene in GENES:
    logging.info(f'Plotting signal near {gene} TSS')
    gene_tss = tss[tss.gene.str.upper()==gene.upper()]
    if len(gene_tss) == 0:
        logging.warning(f'Skipping gene {gene} (not found in TSS file)')
        continue
    assert(gene_tss.chrom.nunique() == 1)
    chrom = gene_tss.chrom.unique()[0]
    start = int(gene_tss.start.min() - FLANK)
    end = int(gene_tss.end.max() + FLANK)
    nearby_tss = tss[(tss.end>=start) & (tss.end<=end)]
    other_tss = nearby_tss[nearby_tss.gene.str.upper()!=gene.upper()]

    fig, axs = plt.subplots(nrows=len(BIGWIGS)+1, ncols=1, figsize=(10, 2*len(BIGWIGS)+1), gridspec_kw={'hspace':0, 'wspace':0, 'height_ratios': [3 for i in BIGWIGS] + [1]})
    
    axs[0].set_title(f'Signal near {gene} TSS ({chrom}:{start}-{end})')

    for ax, i in zip(axs.flatten(), BIGWIGS):
        properties_dict = {'file': i, 'height': 3, 'color': 'black'}
        bw = pygtk.BigWigTrack(properties_dict)
        bw.plot(ax, chrom, start, end)
        ax.set_xlim(left=start, right=end)
        ax.xaxis.set_ticklabels([])
        ax.set_xticks([], [])
        ax.set_ylim(bottom=0)
        ax.legend(title=os.path.basename(i))
    
    max_ylim = max([ax.get_ylim()[1] for ax in axs[:-1]])
    for ax in axs[:-1]:
        ax.set_ylim(0, max_ylim)
        
    # add TSS locations
    ax = axs[-1]
    ax.set_xlim(left=start, right=end)
    for count, i in enumerate(other_tss.end, 1):
        if count == 1:
            ax.axvline(x=i, ymin=0.5, ymax=1, color='black', label='Other gene TSS', ls='--')
        else:
            ax.axvline(x=i, ymin=0.5, ymax=1, color='black', ls='--')
    for count, i in enumerate(gene_tss.end, 1):
        if count == 1:
            ax.axvline(x=i, ymin=0.5, ymax=1, color='red', label=f'{gene} TSS', ls='--')
        else:
            ax.axvline(x=i, ymin=0.5, ymax=1, color='red', ls='--')
    ax.set_yticks([], [])
    ax.legend(loc='right')
    fig.tight_layout()
    fig.savefig(f'{gene}.png', dpi=300)
    fig.clf()

logging.info('Done.')