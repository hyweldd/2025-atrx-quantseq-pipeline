#!/usr/bin/env python

"""
Script to trim fastq.gz file from the start of the first run of 10 A bases.
Written by Hywel Dunn-Davies, February 2025.
Adapted from https://biopython.org/docs/latest/Tutorial/chapter_cookbook.html
"""

import gzip
import argparse
import sys
import os

from collections import defaultdict
from pprint import pprint
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input filepath")
parser.add_argument("-o", "--output", help="Output filepath")
args = parser.parse_args()


def trim_adaptors(records, adaptor):
    """Trims perfect adaptor sequences.

    This is a generator function, the records argument should
    be a list or iterator returning SeqRecord objects.
    """
    len_adaptor = len(adaptor)  # cache this for later
    for record in records:
        index = record.seq.find(adaptor)
        if index == -1: # adaptor not found, so won't trim
            yield record
        else: # trim off the adaptor
            new_record = record[:index]
            yield new_record


count_dict = defaultdict(int)
count = 0

output_fp = args.output
if os.path.exists(output_fp):
    sys.exit(1)

with gzip.open(args.input, "rt") as input_fh:
    with gzip.open(output_fp, "wt") as output_fh:
        original_reads = SeqIO.parse(input_fh, "fastq")
        trimmed_reads = trim_adaptors(original_reads, "AAAAAAAAAA")
        for trimmed_read in trimmed_reads:
            if len(trimmed_read.seq)>0:
                count_dict[len(trimmed_read)] += 1
                count += 1
                output_fh.write(trimmed_read.format("fastq"))

print('_________________')
print(output_fp)
pprint(count_dict)
print("Saved %i reads" % count)
print('_________________')
