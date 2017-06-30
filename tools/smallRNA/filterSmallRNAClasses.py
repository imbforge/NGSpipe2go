#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import (absolute_import, division,
                        print_function)
import argparse
import pysam
# import sys

usage = '''
   Filter classes of small RNAs from library.
   Classes are defined by read sequence length
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
        help='Path to alignments (bam) to be filtered. The bam file needs to be indexed.'
    )

    parser.add_argument(
        '-o', '--outputBam',
        required=True,
        type=str,
        help='Output Bam file to store filtered results.'
    )

    parser.add_argument(
        '-c', '--RNAclass',
        required=True,
        type=str,
        help='Class of small RNAs to be output. Options: 21U, 22G or 26G'
    )

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = getArgs()
    in_file = args.inputBam
    out_file = args.outputBam
    if args.RNAclass == "21U":
        length = 21
        nucleotide_for = "T"
        nucleotide_rev = "A"
    elif args.RNAclass == "22G":
        length = 22
        nucleotide_for = "G"
        nucleotide_rev = "C"
    elif args.RNAclass == "26G":
        length = 26
        nucleotide_for = "G"
        nucleotide_rev = "C"
    else:
        print("Please specify the class of small RNAs")

    piRNA_reads = 0
    other_reads = 0

    inbam = pysam.Samfile(in_file, "rb")
    outbam = pysam.Samfile(out_file, 'wb', template=inbam)

    for read in inbam.fetch():
        if (read.is_reverse is True and
                read.seq.endswith(nucleotide_rev) and
                read.qlen == length):
            piRNA_reads = piRNA_reads + 1
            outbam.write(read)
        if (read.is_reverse is False and
                read.seq.startswith(nucleotide_for) and
                read.qlen == length):
            piRNA_reads = piRNA_reads + 1
            outbam.write(read)
        else:
            other_reads = other_reads + 1

    outbam.close()
    inbam.close()

    pysam.index(out_file)

    print("Number of %s reads: %d" % (args.RNAclass, piRNA_reads))
    print("Number of reads filtered out: %d" % other_reads)
    outbam.close()
