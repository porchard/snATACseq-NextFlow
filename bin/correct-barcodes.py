#!/usr/bin/env python

import os
import gzip
import sys
import logging
from multiprocessing import Pool
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import functools 
import argparse 

parser = argparse.ArgumentParser(description='Given a set of barcode fastq files and the 10X barcode whitelist, output barcode corrections (read name, old barcode, new barcode).', add_help = True)
parser.add_argument('--threads', type = int, default = 1, help = 'Number of threads to use.')
parser.add_argument('whitelist', type = str, default = '', help = '10X whitelist file')
parser.add_argument('fastq', type = str, nargs = '+', default = '', help = 'Fastq files')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

PREFIX_SIZE = 10

def hamming_distance(x, y):
    hamming_dist = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            hamming_dist += 1
    return hamming_dist


def hamming_distance_within_n(n, x, y, start_from = 0):
    hamming_dist = 0
    for i in range(start_from, len(x)):
        if x[i] != y[i]:
            hamming_dist += 1
            if hamming_dist > n:
                return False
    return True


# read in the barcode fastqs...
logging.info('Reading fastq files to get priors...')
counts = dict()

for fastq in args.fastq:
    logging.info(f'Processing fastq file {fastq}')
    read_count = 0
    with gzip.open(fastq, 'rt') as f:
        for read_name, uncorrected_barcode, phred in FastqGeneralIterator(f):
            read_count += 1
            if read_count % 1000000 == 0:
                logging.info(f'Processed {read_count} reads')
            if uncorrected_barcode not in counts:
                counts[uncorrected_barcode] = 0
            counts[uncorrected_barcode] += 1

logging.info('Finished fetching initial barcode counts.')


# read in the whitelist
logging.info('Reading whitelist')
good_barcodes = set()
with open(args.whitelist, 'r') as f:
    for line in f:
        b = line.rstrip()
        good_barcodes.add(b)
        if b not in counts:
            counts[b] = 0
        counts[b] += 1 # add pseudocount
        

logging.info('Finished reading whitelist')

# I believe cellranger adds a pseudocount
# to unobserved whitelisted reads so that barcodes
# can still be corrected to unobserved whitelisted reads.
# I will not to this though.
bad_barcodes = list(set([b for b in counts if b not in good_barcodes]))
logging.info('Found {} unique whitelisted barcodes'.format(len(counts) - len(bad_barcodes)))
logging.info('Found {} unique non-whitelisted barcodes'.format(len(bad_barcodes)))

logging.info('Pre-calculating hamming distances between non-whitelisted barcodes and observed whitelisted barcodes')

logging.debug('Building whitelist tree')
whitelist_tree = dict() # first PREFIX_SIZE bp --> whitelist
for i in good_barcodes:
    prefix = i[0:PREFIX_SIZE]
    if prefix not in whitelist_tree:
        whitelist_tree[prefix] = []
    whitelist_tree[prefix].append(i)

logging.debug('Building barcode tree')
barcode_tree = dict()
prefixes = set([b[0:PREFIX_SIZE] for b in bad_barcodes])

def add_prefix(prefix):
    tree = []
    for p in whitelist_tree:
        if hamming_distance_within_n(2, prefix, p):
            tree.append([p, hamming_distance(prefix, p)])
    return (prefix, tree)

# add in the prefix
x = None
with Pool(args.threads) as p:
    x = p.map(add_prefix, list(prefixes))
for key, val in x:
    barcode_tree[key] = val

logging.debug('Finished building barcode tree')

# find the true hamming distance ones for each barcode
logging.info('Finding the hamming distances.')
possible_fixes = dict() # bad barcode --> whitelisted ones within hamming distance of 2

def add_possible_fix(barcode):
    fixes = []
    for prefix, dist in barcode_tree[barcode[0:PREFIX_SIZE]]:
        for fix in whitelist_tree[prefix]:
            if hamming_distance_within_n(2 - dist, barcode, fix, start_from = PREFIX_SIZE):
                fixes.append(fix)
    return (barcode, fixes)

with Pool(args.threads) as p:
    x = p.map(add_possible_fix, bad_barcodes)
for key, val in x:
    possible_fixes[key] = val

logging.info('Finished adding possible fixes')


def prob_of_error(uncorrected, potential_correction, phred):
    # find the error
    total_p = 1
    for i in range(len(uncorrected)):
        if uncorrected[i] != potential_correction[i]:
            q = ord(phred[i]) - 33
            p = 10**(-q/10)
            total_p *= p
    return total_p if total_p != 1 else None

def correct_barcode(uncorrected, potential_fix, phred):
    if len(potential_fix) == 1:
        return potential_fix[0]
    elif len(potential_fix) == 0:
        return None
    else:
        priors = [counts[barcode] for barcode in potential_fix]
        prob_of_errors = [prob_of_error(uncorrected, barcode, phred) for barcode in potential_fix]
        relative_magnitude = [priors[i] * prob_of_errors[i] for i in range(len(priors))]
        if sum(relative_magnitude) == 0:
            return None
        relative_magnitude = [i/sum(relative_magnitude) for i in relative_magnitude]
        if max(relative_magnitude) >= 0.975:
            for i in range(len(potential_fix)):
                if relative_magnitude[i] >= 0.975:
                    return potential_fix[i]
        else:
            return None


logging.info('Reading through fastq files again to generate fixes...')

fixes = dict()
for fastq in args.fastq:
    logging.info(f'Processing fastq file {fastq}')
    read_count = 0
    with gzip.open(fastq, 'rt') as f:
        for read_name, uncorrected_barcode, phred in FastqGeneralIterator(f):
            read_count += 1
            if read_count % 1000000 == 0:
                logging.info(f'Processed {read_count} reads')
            if uncorrected_barcode not in good_barcodes:
                read_name = read_name.split(' ')[0]
                key = uncorrected_barcode + ':' + phred
                corrected_barcode = correct_barcode(uncorrected_barcode, possible_fixes[uncorrected_barcode], phred) if key not in fixes else fixes[key]
                fixes[key] = corrected_barcode
                print('{}\t{}\t{}'.format(read_name, uncorrected_barcode, corrected_barcode))

logging.info('Done.') 
