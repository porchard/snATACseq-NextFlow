#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import re
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('bam_in', help = 'Bam file to adjust tags in')
parser.add_argument('bam_out', help = 'Bam file to write')
args = parser.parse_args()


with pysam.AlignmentFile(args.bam_in, 'rb') as f:
    with pysam.AlignmentFile(args.bam_out, 'wb', template = f) as f_new:
        for read in f.fetch(until_eof=True):
            read.set_tag('RG', read.get_tag('CB'), replace = False)
            f_new.write(read)

