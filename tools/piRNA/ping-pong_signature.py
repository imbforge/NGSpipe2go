#!/usr/bin/env python
# encoding: utf-8


usage = '''
    Calculates the overlap between reads (bam), assumed to be piRNAs, using pybedtools. Something like this:

    5'-----------3' primary
        3'-------------5' secondary

    If an list of genomic features is provided (in bed format), only reads overlapping with features will be analysed. Otherwise all reads are analysed.
    Only reads with a certain size are analysed (default > 21 & < 34)
    When all reads are used, those in the "+" strand are defined as "primary". If an a list of genomic regions is provided, then one needs to define which reads are "primary". By default, those which are in the same strand as the interval will be defined as "primary".
    Intersection is by 1bp and strand-ware, i.e., reads must be in opposite strands to be counted.
'''

import argparse
import os
import pybedtools
from pybedtools import BedTool
import sys


def getArgs():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )

    parser.add_argument(
        '-b', '--bam',
        required=True,
        type=str,
        help='Path to alignments (bam) to look for motif.'
    )

    parser.add_argument(
        '-i', '--intervals',
        required=False,
        type=str,
        help='Path to genomic intervals (bed), usually containing the genomic locations of repeat elements. This will be used to determine sense or antisense of read mapping'
    )

    parser.add_argument(
        '--minsize',
        required=False,
        type=int,
        default=21,
        help='Smallest size of reads to be analysed. Default is 21 nucleotides'
    )

    parser.add_argument(
        '--maxsize',
        required=False,
        type=int,
        default=34,
        help='Largest size of reads to be analysed. Default is 34 nucleotides.'
    )

    parser.add_argument(
        '--primary',
        required=False,
        type=str,
        default="sense",
        help='Which reads are to be considered "primary". A choice of "sense" or "antisense" relative to the provided genomic regions (--intervals). When no genomic regions are provided "sense/antisense" stand for "+/-" strands. Default: "sense" (+)'
    )

    parser.add_argument(
        '--maxres',
        required=False,
        type=int,
        default=30,
        help='Maximal overlapping distance to be reported. Serves only to avoid extra data. Not used yet. Default: 30'
    )

    parser.add_argument(
        '-o', '--outFolder',
        required=False,
        type=str,
        default="./",
        help='Folder to put the results. Default: current folder (./)'
    )

    args = parser.parse_args()
    return args


def createAndChangeDir(dir_name):
    '''
    Creates a directory from a string (dir_name) and changes to that directory.
    '''
    try:
        os.makedirs(dir_name)
    except OSError:
        if not os.path.isdir(dir_name):
            raise
    os.chdir(dir_name)


def filterReadsByLength(inbam, minlength, maxlength):
    '''
    Takes a bam file and selects intervals that are within the defined lengths.
    Input: bam file and min/max lengths
    Output: bedTool
    '''
    # convert bam to bed
    intervals = BedTool(inbam).bam_to_bed()
    filt = intervals.filter(lambda x: len(x) > minlength and len(x) < maxlength).saveas()
    # print filt
    return filt


def separateStrands(reads):
    '''
    Splits bedTools intervals in plus and minus strands.
    Input: bedTools object
    Output: list of tuples with a name (str) and the reads for each strand (bedTool)
    '''
    # # convert bam to bed
    # piRNA = BedTool(inbam).bam_to_bed()

    antisense = reads.filter(lambda b: b.strand == '-').saveas()
    sense = reads.filter(lambda b: b.strand == '+').saveas()
    strands = [
        ('sense', sense),
        ('antisense', antisense)
        ]
    return strands


def intersectReadsWithFeatures(reads, inbed):
    '''
    Intersects reads with genomic features, e.g. Transposable elements, and returns separately reads that map sense and antisense to the features.
    Input: a bedtools object and path to a bed file
    Output: list of tuples with a name (str) and the reads for sense and antisense piRNAs (bedTool)
    '''

    ## create bedtool for genomic features
    bed = BedTool(inbed)
    ## separate sense and antisense reads
    sense = reads.intersect(bed, s=True)
    antisense = reads.intersect(bed, S=True)
    piRNAs = [
        ('sense', sense),
        ('antisense', antisense)
    ]
    return piRNAs


def intersectIntervals(sense, antisense):
    '''
    Intersects two list of genomic intervals, and returns intersecting intervals.
    Input: bedtools intervals
    Output:
    '''
    pairs = sense.intersect(antisense, wo=True, S=True)
    print len(pairs)
    return pairs


def overlapDistance(feature):
    '''
    Calculates the 5'-5' difference between pair of overlaps.
    Input: bedTool consisting of both a and b from an intersection. Assumes bed6
    '''
    if feature[5] == "+":
        pri_p5 = feature[1]
        sec_p5 = feature[8]
        dist = int(sec_p5) - int(pri_p5)
    elif feature[5] == "-":
        pri_p5 = feature[2]
        sec_p5 = feature[7]
        dist = int(pri_p5) - int(sec_p5)
    return dist


def countFrequency(lst):
    '''
    Calculates the frequency of each 5'-5' interval
    '''
    d = {}
    for item in lst:
        if item in d:
            d[item] = d.get(item)+1
        else:
            d[item] = 1

    # for k, v in d.items():
    #     print(str(k)+'\t'+str(v))
    return d


def writeResults(d, filename):
    with open(filename, 'w') as f:
        for k, v in d.items():
            line = '{}\t{}\n'.format(k, v)
            f.write(line)


def computePairOverlap(strands, primary, maxres):
    '''
    Wrapper function that takes sense and antisense strands and calculates the overlap
    Input: Tupple with sense and antisense strands.
    output: a file with Position in the first col and frequency in the second
    '''
    # print("Primary is set to {}".format(primary))
    strds = dict(strands)
    if primary == "sense":
        print 'primary is sense/+'
        pairs = intersectIntervals(
            strds["sense"],
            strds["antisense"]
            )
    elif primary == "antisense":
        print 'primary is antisense/-'
        pairs = intersectIntervals(
            strds["antisense"],
            strds["sense"]
            )

    ## "loops" over each feature and calculates distance
    pp = map(overlapDistance, pairs)
    freq = countFrequency(pp)
    writeResults(freq, "pp_freq.txt")


if __name__ == '__main__':
    args = getArgs()
    print args
    inbam = args.bam
    inbed = args.intervals
    minlength = args.minsize
    maxlength = args.maxsize
    primary = args.primary
    maxres = args.maxres
    out_folder = args.outFolder

    ## filter reads too short or long:
    reads = filterReadsByLength(inbam, minlength, maxlength)
    # print(reads)
    ## define sense and antisense
    if inbed is None:
        print 'No genomic regions have been provided. Considering sense and antisense based on read strand.'
        strands = separateStrands(reads)
        # print strands[1][1]
    else:
        print 'Intersecting with genomic features to define sense and antisense piRNAs.'
        strands = intersectReadsWithFeatures(reads, inbed)

    createAndChangeDir(out_folder)
    computePairOverlap(strands, primary, maxres)
    # os.system("Rscript ~/imb-kettinggr/adomingues/projects/ping-pong/scripts/plotPP.R")
    script_dirname = os.path.dirname(os.path.realpath(sys.argv[0]))
    plot_script_path = script_dirname + '/plotPP.R'
    os.system("Rscript " + plot_script_path)
