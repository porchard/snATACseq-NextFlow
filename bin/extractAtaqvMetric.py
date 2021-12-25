#!/usr/bin/env python

import argparse
import json
import gzip
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--metrics', nargs = '*', default = None, help = 'Names of metrics to extract. If not given, print all numeric metrics.')
parser.add_argument('--files', nargs = '+', required=True, default = [], help = 'Names of ataqv metric files from which to extract metrics. These files may be output by ataqv itself, or be within the data/ directory of an ataqv app.')
args = parser.parse_args()

def median_from_dict(d):
    """
    Given a dictionary of value counts (d[value] = count), calculate the median value.

    d: dictionary, where values are non-negative integers.
    """
    if not isinstance(d, dict):
        raise TypeError('d must be a dict')
    if not len(d) > 0:
        raise ValueError('Cannot find median of an empty dict')
    for i in d.values():
        if not isinstance(i, (int, np.integer)):
            raise TypeError('All values in dict d must be integers')
        if i < 0:
            raise ValueError('All values in dict d must be >= 0')
    count = sum(d.values())
    middle = [int(count / 2), int(count / 2) + 1] if count % 2 == 0 else [int(count / 2) + 1]
    middle_values = []
    index = 0
    for key, val in sorted(d.items()):
        for i in middle:
            if i > index and i <= (index + val):
                middle_values.append(key)
        index += val
        if len(middle_values) == len(middle):
            break
    return np.mean(middle_values)


def _add_metrics(d):
    d['percent_mitochondrial'] = 100 * float(d['total_mitochondrial_reads']) / d['total_reads'] if d['total_reads'] > 0 else None
    d['percent_mitochondrial_duplicate'] = 100 * float(d['duplicate_mitochondrial_reads']) / d['total_mitochondrial_reads'] if d['total_mitochondrial_reads'] > 0 else None
    d['percent_autosomal_duplicate'] = 100 * float(d['duplicate_autosomal_reads']) / d['total_autosomal_reads'] if d['total_autosomal_reads'] > 0 else None
    d['percent_hqaa'] = 100.0 * float(d['hqaa']) / d['total_reads']
    d['percent_properly_paired_and_mapped'] = 100.0 * float(d['properly_paired_and_mapped_reads']) / d['total_reads']
    d['percent_secondary'] = 100.0 * float(d['secondary_reads']) / d['total_reads']
    d['percent_supplementary'] = 100.0 * float(d['supplementary_reads']) / d['total_reads']
    d['percent_duplicate'] = 100.0 * float(d['duplicate_reads']) / d['total_reads']
    d['percent_unmapped'] = 100.0 * float(d['unmapped_reads']) / d['total_reads']
    d['percent_unmapped_mate'] = 100.0 * float(d['unmapped_mate_reads']) / d['total_reads']
    d['percent_qcfailed'] = 100.0 * float(d['qcfailed_reads']) / d['total_reads']
    d['percent_unpaired'] = 100.0 * float(d['unpaired_reads']) / d['total_reads']
    d['percent_mapq_0'] = 100.0 * float(d['reads_mapped_with_zero_quality']) / d['total_reads']
    d['percent_rf'] = 100.0 * float(d['rf_reads']) / d['total_reads']
    d['percent_ff'] = 100.0 * float(d['ff_reads']) / d['total_reads']
    d['percent_rr'] = 100.0 * float(d['rr_reads']) / d['total_reads']
    d['percent_autosomal'] = 100.0 * float(d['total_autosomal_reads']) / d['total_reads'] if d['total_reads'] > 0 else None
    d['percent_mate_separate_chromosome'] = 100.0 * float(d['reads_with_mate_mapped_to_different_reference']) / d['total_reads']
    d['percent_mate_too_distant'] = 100.0 * float(d['reads_with_mate_too_distant']) / d['total_reads']
    d['percent_improperly_paired'] = 100.0 * float(d['reads_mapped_and_paired_but_improperly']) / d['total_reads']
    d['median_fragment_length'] = median_from_dict({i[0]: i[1] for i in d['fragment_length_counts']})
    return d


def read_ataqv_json(metrics_file):
    """
    Load ataqv metrics from an ataqv JSON file. Also adds some metrics that may not otherwise be present (e.g., median_fragment_length).

    metrics_file: path to JSON file.

    Return: list of dictionaries.
    """
    with gzip.open(metrics_file, 'rb') as f:
        line = f.read()
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        l = json.loads(line)
        if 'metrics' in l[0]:
            l = [d['metrics'] for d in l]
        l = [_add_metrics(d) for d in l]
        return l


for ataqv_data_file in args.files:
    for d in read_ataqv_json(ataqv_data_file):
        if args.metrics is None:
            for metric, value in d.items():
                if isinstance(value, (float, int)):
                    print('{}\t{}\t{}\t{}'.format(ataqv_data_file, d['name'], metric, value))
        else:
            for metric in args.metrics:
                if not metric in d:
                    raise ValueError('Metric {} not present'.format(metric))
                print('{}\t{}\t{}\t{}'.format(ataqv_data_file, d['name'], metric, d[metric]))