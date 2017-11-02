#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (absolute_import, division,
                        print_function)
import argparse
import pysam
import sys

usage = '''
   Filter classes of small RNAs from library.
   Classes are defined by read sequence length
   '''


def getArgs():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-i', '--inputBam',
        required=True,
        type=str,
        help='Path to alignments (bam) to be filtered. The bam file needs to be indexed.'
    )

    parser.add_argument(
        '-o', '--outputBam',
        required=True,
        type=str,
        help='Output Bam file to store filtered results. Use "stdout" to write output to standard out'
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
        '-n', '--nuc',
        required=False,
        type=str,
        help='Alignments with nuc at the first position, or the reverse complement at the the last position if read is mapped in the reverse strand, will be kept.'
    )

    args = parser.parse_args()

    if not (args.min or args.max or args.nuc):
        parser.error(
            'Filter option not set. Add at least one of the following filtering options:\n --min, --max or --nuc'
            )
    return args


def eprint(*args, **kwargs):
    """
    helper function to print to stderr
    source:
    https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


if __name__ == '__main__':
    args = getArgs()
    in_file = args.inputBam
    out_file = args.outputBam

    if args.nuc:
        if args.nuc == "T":
            nucleotide_for = "T"
            nucleotide_rev = "A"
        elif args.nuc == "A":
            nucleotide_for = "A"
            nucleotide_rev = "T"
        elif args.nuc == "G":
            nucleotide_for = "G"
            nucleotide_rev = "C"
        elif args.nuc == "C":
            nucleotide_for = "C"
            nucleotide_rev = "G"
        else:
            eprint("Not a valid nucleotide. Use one of: A, C, G, T")


    inbam = pysam.AlignmentFile(in_file, "rb")

    if out_file == "stdout":
        outbam = pysam.AlignmentFile("-", "wb", template=inbam)
    else:
        outbam = pysam.AlignmentFile(out_file, 'wb', template=inbam)

    reads_total = 0
    reads_kept = 0
    for read in inbam.fetch():
        reads_total = reads_total + 1
        if args.min:
            if read.qlen < args.min:
                continue
        if args.max:
            if read.qlen > args.max:
                continue
        if args.nuc:
            if (
                read.is_reverse is True and
                read.seq.endswith(nucleotide_rev) is False):
                continue
            if (
                read.is_reverse is False and
                read.seq.startswith(nucleotide_for) is False):
                continue
        reads_kept = reads_kept + 1
        outbam.write(read)
    inbam.close()

    if out_file != "stdout":
        outbam.close()
        pysam.index(out_file)

    eprint("Total number of reads in input: %d" % reads_total)
    eprint("Number of reads after filtering: %d" % (reads_kept))
