#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (absolute_import, division,
                        print_function)
import argparse
import pysam
import string
import sys

usage = '''
   Filter classes of small RNAs from library.
   Classes are defined by read sequence length and/or a nucleotide in a particular sequence position.
   NOTE: assumes bowtie mapping, and thus sequence alignments in the reverse strand will be reverse complemented before testing for nucleotide position. The output bam will keep the same sequence alignment as the original bam.
   '''


def getArgs():
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-i', '--inputBam',
        required=True,
        type=str,
        help='Path to alignments (bam) to be filtered. The bam file needs to be indexed. Use "stdin" to read input from standard input (stream). If stdin the file is assumed to be sam and the header is included. Example: samtools view -h data/ah.bam | filterReads/filterSmallRNAclasses.py -i stdin -m 21 -M21 -n "A" -p last -o stdout'
    )

    parser.add_argument(
        '-o', '--outputBam',
        required=True,
        type=str,
        help='Output Bam file to store filtered results. Use "stdout" to write output to standard out (stream).'
    )

    parser.add_argument(
        '-m', '--min',
        required=False,
        type=int,
        help='Reads with length >= m will be keep.'
    )

    parser.add_argument(
        '-M', '--max',
        required=False,
        type=int,
        help='Reads with length <= M will be keep.'
    )

    parser.add_argument(
        '-l', '--length',
        required=False,
        type=int,
        help='The exact length (l) of reads to be kept. Not implemented'
    )

    parser.add_argument(
        '-n', '--nucleotide',
        required=False,
        type=str,
        choices=('T', 'G', 'A', 'C'),
        help='Filter by nucleotide'
    )

    parser.add_argument(
        '-p', '--position',
        required=False,
        type=str,
        choices=('first', 'last'),
        help='Should the nucleotide be in the first or last sequence position (5\' or 3\' end). Here sequence is considered to be original sequenced DNA, so if the is mapped in the reverse strand the alignment sequence will be reverse completed.'
    )

    args = parser.parse_args()

    if not (args.min or args.max or args.nucleotide):
        parser.error(
            'Filter option not set. Add at least one of the following filtering options:\n --min, --max or --nucleotide'
            )

    if args.nucleotide and not args.position:
        parser.error("Filtering by nucleotide requires the position argument.")

    return args


def eprint(*args, **kwargs):
    """
    helper function to print to stderr
    source:
    https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


def reverseComplement(seq):
    # credit to: https://bioinformatics.stackexchange.com/a/3585/108
    tab = str.maketrans("ACTG", "TGAC")
    return seq.translate(tab)[::-1]


if __name__ == '__main__':
    args = getArgs()
    in_file = args.inputBam
    out_file = args.outputBam
    min_length = args.min
    max_length = args.max
    nucleotide = args.nucleotide
    position = args.position


    if in_file == "stdin":
        inbam = pysam.AlignmentFile("-", "r")
    else:
        inbam = pysam.AlignmentFile(in_file, "rb")

    if out_file == "stdout":
        outbam = pysam.AlignmentFile("-", "wb", template=inbam)
    else:
        outbam = pysam.AlignmentFile(out_file, 'wb', template=inbam)

    reads_total = 0
    reads_kept = 0

    for read in inbam.fetch():
        reads_total = reads_total + 1
        if min_length:
            if read.qlen < min_length:
                continue
        if max_length:
            if read.qlen > max_length:
                continue

        if args.nucleotide:
            if read.is_reverse is True:
                sequence = reverseComplement(read.seq)
            else:
                sequence = read.seq

            if (
                position == "first" and
                sequence.startswith(nucleotide) is False
            ):
                continue

            if (
                position == "last" and
                sequence.endswith(nucleotide) is False
            ):
                continue

        reads_kept = reads_kept + 1
        outbam.write(read)
    inbam.close()

    if out_file != "stdout":
        outbam.close()
        pysam.index(out_file)

    eprint("Total number of reads in input: %d" % reads_total)
    eprint("Number of reads after filtering: %d" % (reads_kept))
