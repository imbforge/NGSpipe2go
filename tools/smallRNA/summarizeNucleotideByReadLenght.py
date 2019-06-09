#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (absolute_import, division,
                        print_function)
import argparse
import pysam
from collections import Counter
# import sys


usage = '''
   Counts the frequency of each nucleotide per read length.
   '''


def getArgs():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-i', '--inputBam',
        required=True,
        type=str,
        help='Path to alignments (bam) to be summarized. The index needs to be present. Use "stdin" to read input from standard input (stream). If stdin the file is assumed to be sam and the header is included.'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        help='Path to file to store summary.'
    )

    args = parser.parse_args()
    return args


def getComplement(base):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complement = basecomplement[base.upper()]
    return complement


def writeResults(res, file_out):
    with open(file_out, 'w') as f:
        for k, v in res:
            st = "{}\t{}\t{}".format(k[0], k[1], v)
            f.write(st + '\n')


if __name__ == '__main__':
    args = getArgs()
    in_file = args.inputBam
    out_file = args.output

    if in_file == "stdin":
        inbam = pysam.AlignmentFile("-", "r")
    else:
        inbam = pysam.AlignmentFile(in_file, "rb")

    list_of_pairs = []
    for read in inbam.fetch():
        if (read.is_reverse is True):
            nucleotide = getComplement(read.seq[-1])
        if read.is_reverse is False:
            nucleotide = read.seq[0]
        read_len = read.query_length

        list_of_pairs.append((read_len, nucleotide))
    inbam.close()

    stats = Counter(list_of_pairs).most_common()

    writeResults(stats, out_file)
