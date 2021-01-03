#!/usr/bin/env python
# encoding: utf-8

usage = '''

Takes alignments, corresponding to piRNAs, and counts how often each nucleotide occurs in the boundaries (5' and 3') of the sequence. By default the it counts nucleotide frequency at positions -20 to +20. It also outputs results for sense and antisense piRNAs to a provided list of genomic features. A plot in pdf and png format is generated along with the files containing the counts.

Author: António Domingues
amjdomingues [at] gmail.com
'''


import pybedtools
from pybedtools import BedTool
from pybedtools.featurefuncs import three_prime, five_prime, greater_than
import csv
import os
import sys
import argparse
import pandas as pd
import pysam


def getArgs():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-f', '--fasta',
        required=True,
        type=str,
        help='Path to fasta to retrieve sequences.'
    )

    parser.add_argument(
        '-b', '--bam',
        required=True,
        type=str,
        help='Path to alignments (bam) to look for motif.'
    )

    parser.add_argument(
        '-i', '--intervals',
        required=True,
        type=str,
        help='Path to genomic intervals (bed), usually containing the genomic locations of repeat elements. This will be used to determine sense or antisense of read mapping'
    )

    parser.add_argument(
        '-u', '--upstream',
        required=False,
        type=int,
        default=20,
        help='Number of nucleotides upstream of the read (genomic sequence). Default is 20'
    )

    parser.add_argument(
        '-d', '--downstream',
        required=False,
        type=int,
        default=20,
        help='Number of nucleotides downstream of the read (piRNA sequence). Default is 20'
    )

    parser.add_argument(
        '-g', '--genome',
        required=False,
        type=str,
        help='The genome to retrieve chromosome lengths, e.g., hg19, mm10, danRer7...'
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


def timeStamp():
    '''
    Returns current system time in the format "YYYY-MM-DD HH:MM:SS"
    Source: http://stackoverflow.com/a/13891070/1274242
    '''
    import time
    import datetime
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return(st)


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


def get_chrom_lengths(path_to_bam):
    '''
    Uses pysam to retrieve chromosome sizes from bam.
    Useful helper to use with some pybedtools functions (e.g. coverage), when a bam was mapped with custom genome not available in UCSC.
    Input: path to bam file (should be indexed)
    Output: dictionary.
    Example output:
    {'chr4': (0, 1351857), 'chr3L': (0, 24543557), 'chr2L': (0, 23011544), '*': (0, 0), 'chrX': (0, 22422827), 'chr2R': (0, 21146708), 'chr3R': (0, 27905053)}
    '''
    idx = pysam.idxstats(path_to_bam).splitlines()
    chromsizes = {}
    for element in idx:
        stats = element.split("\t")
        chromsizes[stats[0]] = (0, int(stats[1]))
    return chromsizes


def countNucleotidePerPosition(sequences):
    '''
    Takes a list of strings and determines the nucleotide occurrence per position.
    Returns a panda DataFrame where each column is a nucleotide and each row a position.

    source: http://stackoverflow.com/a/21103385/1274242
    '''
    print('Counting nucleotides', timeStamp())
    df = pd.DataFrame([list(s) for s in sequences])
    counts = df.apply(pd.value_counts).transpose()
    return(counts)


def parseSequence(entries):
    '''
    Extracts the nucleotide sequence from the poorly formated output of BedTool.sequence.
    Input example:
    >chrDESTTol2CG2:2717-2727\nCAAAATCCAC\n>chrDESTTol2CG2:3713-3723\nGTCGATGCCC\n
    '''
    seqs = []
    for entry in entries.split('\n'):
        entry = entry.strip('\n').strip(' ')
        if '>' in entry:
            pass
        else:
            seqs.append(entry.upper())
    seqs = filter(None, seqs)
    return seqs


def getSequences(coordinates, fasta):
    '''
    Retrieves the sequence corresponding to a coordinate.
    Input: BedTools object
    Output:
    >chrDESTTol2CG2:2717-2727\nCAAAATCCAC\n>chrDESTTol2CG2:3713-3723\nGTCGATGCCC\n
    '''
    sequence = coordinates.sequence(fi=fasta, s=True)
    return open(sequence.seqfn).read()


def getSequencesFrom5prime(coordinates, fasta, upstream, downstream, chromsizes):
    '''
    Retrieves the sequences limiting the 5' end of a genomic location.
    Input: BedTools object
    Output: a list of strings with defined size surrounding a genomic location
    '''
    print('Fetching 5\' sequences', timeStamp())
    seq_len = upstream + downstream - 1
    start = coordinates.each(five_prime, upstream=upstream, downstream=downstream,
                             genome=chromsizes).filter(greater_than, seq_len).saveas()
    clean_seq = parseSequence(getSequences(start, fasta))
    return clean_seq


def getSequencesFrom3prime(coordinates, fasta, upstream, downstream, chromsizes):
    '''
    Retrieves the sequences limiting the 3' end of a genomic location.
    Input: BedTools object
    Output: a list of strings with defined size surrounding a genomic location
    '''
    print('Fetching 3\' sequences', timeStamp())
    seq_len = upstream + downstream - 1
    start = coordinates.each(three_prime, upstream=upstream, downstream=downstream,
                             genome=chromsizes).filter(greater_than, seq_len).saveas()
    clean_seq = parseSequence(getSequences(start, fasta))
    return clean_seq


def convertSGVImages(image, res=300):
    '''
    Convert an svg image into a high-resolution png
    Input: path to image
    '''
    image_out = image.replace('.svg', '.png')
    command = 'convert -density ' + str(res) + ' ' + image + ' ' + image_out
    print('converting logo to PNG:', timeStamp())
    os.system(command)


# def createMotif(sequences, fname):
#     '''
#     Creates and saves a motif (logo) from a list of input sequences.
#     Input: list of strings with the sequences; file path to save the logo.
#     Output: an image file with the logo.
#     '''
#     print('Generating motif', timeStamp())
#     try:
#         os.makedirs('figure')
#     except OSError:
#         if not os.path.isdir('figure'):
#          raise
#     from Bio.Seq import Seq
#     from Bio import motifs
#     from Bio.Alphabet import IUPAC
#     import urllib2

#     # m = motif.motif(alphabet=IUPAC.unambiguous_dna) # initialize motif
#     instances = []
#     for sequence in sequences:
#         if len(sequence) < 40:
#             print(sequence)
#         instances.append(Seq(sequence, alphabet=IUPAC.ambiguous_dna))
#     m = motifs.create(instances)
#     flogo = 'figure/' + fname
#     while True:
#         # source: http://stackoverflow.com/a/9986206/1274242
#         try:
#             m.weblogo(flogo, format='SVG')
#             break
#         except urllib2.HTTPError, detail:
#             if detail.errno == 500:
#                 time.sleep(5)
#                 continue
#             else:
#                 raise
#     convertSGVImages(flogo)


def intersectBamWithBed(inbam, inbed):
    '''
    Intersects reads with genomic features, Transposable elements, and returns separately reads that map sense and antisense to the features.
    Input: paths to bam and bed file
    Output: list of tuples with a name (str) and the reads for sense and antisense piRNAs (bedTool)
    '''
    # convert bam to bed
    print('Separating sense and antisense piRNAs', timeStamp())
    piRNA = BedTool(inbam).bam_to_bed()
    # create bedtool for genomic features
    bed = BedTool(inbed)   
    # outname = inbam.replace('.bam', '')
    # outsense = outname + "sense.bed"
    # outantisense = outname + "antisense.bed"
    antisense = piRNA.intersect(bed, S=True)
    sense = piRNA.intersect(bed, s=True)
    piRNAs = [
        ('sense', sense),
        ('antisense', antisense)]
    return piRNAs


if __name__ == '__main__':
    args = getArgs()
    print(args)
    fasta = BedTool(args.fasta)
    inbam = args.bam
    inbed = args.intervals
    up = args.upstream
    down = args.downstream
    out_folder = args.outFolder
    if args.genome is not None and not "none":
        chromsizes = pybedtools.chromsizes(args.genome)
    else:
        print('Retrieving custom chromosome lengths', timeStamp())
        chromsizes = get_chrom_lengths(inbam)

    exp_name = os.path.split(inbam)[1].split(
        ".")[0] + os.path.split(inbed)[1].split(".")[0]
    exp_folder = out_folder + '/' + exp_name
    createAndChangeDir(out_folder)
    print('Starting analysis for', inbam, timeStamp())

    strand_piRNAs = intersectBamWithBed(inbam, inbed)
    for direction in strand_piRNAs:
        print('piRNAs in', direction[0], timeStamp())
        print('Fetching 5\' sequences', timeStamp())
        five_out = direction[0] + '.5prime.count.csv'
        five = getSequencesFrom5prime(
            direction[1], fasta, upstream=up, downstream=down, chromsizes=chromsizes)
        
        print('Fetching 3\' sequences', timeStamp())
        three_out = direction[0] + '.3prime.count.csv'
        three = getSequencesFrom3prime(
            direction[1], fasta, upstream=up, downstream=down, chromsizes=chromsizes)
        
        # print 'Counting in parallel ' + timeStamp()
        # pool = multiprocessing.Pool(processes=2)
        # five_count, three_count = pool.map(countNucleotidePerPosition,
        #   (five, three))
        # pool.close()
        # pool.join()
        
        five_count = countNucleotidePerPosition(five)
        index = list(range(-up, 0) + range(1, down + 1))
        five_count['Position'] = index
        five_count.to_csv(five_out, index=False)
        
        three_count = countNucleotidePerPosition(three)
        # reversed_index = [i * -1 for i in list(reversed(range(-up, 0) + range(1, down+1)))]
        three_count['Position'] = list(reversed(index))
        three_count.to_csv(three_out, index=False)
        
        # create motifs:
        # m_out_five = direction[0] + '.5prime.logo.svg'
        # createMotif(five, m_out_five)
        
        # m_out_three = direction[0] + '.3prime.logo.svg'
        # createMotif(three, m_out_three)
    print('Plotting results', timeStamp())
    script_dirname = os.path.dirname(os.path.realpath(sys.argv[0]))
    plot_script_path = script_dirname + '/piRNABaseTerminalBasesPlot.R'
    os.system("Rscript " + plot_script_path)
