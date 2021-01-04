# piRNA pipeline

This pipeline is geared toward zebrafish, but there is no reason to believe that it will not work for other species/taxa which share similar piRNA biogenesis pathtways (_Homo sapiens_, _Mus musculus_, _Bombix morix_, etc). The final results are plots summarizing:

- ping-pong strength
- nucleotide bias, including up and downstream (see Mohn et al), though it could be in a better plot
- profile of targeted RNA species (tranposons)

It will also generate:

- mapped `bam` files with multimapping and unique reads separately
- bigwigs of uniquely mapping reads, with both strands or separating the strands.


# Inputs

- raw reads
- genome sequence (fasta)
- genome annotation (gtf)
- rRNA sequences
- repeatmasker annotations

For examples of how to set-up the files see the [test data folder](testdata/piRNA/README.md)


# Main steps and parameters

- QC
    + adapter trimming
    + removal of reads with low-quality bases
    + removal of PCR duplicates
    + FastQC, performed before and after mapping
    + diagnostics plots
- mapping with bowtie, allowing for one randomly selected alignment for multi-mapping reads, and two mismatches
- ping-pong signal overall ping-pong strength in the library (useful to address if mutants affect piRNA biogenesis)
- [PingPongPro](https://github.com/suhrig/pingpongpro/) to detect in which transposons the ping-pong pathways is active
- nucleotide bias of small RNAs
- counting and summarization of reads mapping to transposable elements to determine if there is a strand  or length bias in particular classes of transposons.
- generate bigwigs 


# Limitations and future improvements

- sense and antisense bigwigs are internally normalized which means that the sum of both signals does not equal the total. Ideally they would be normalized to the total mapped reads.
- it is possible that some double counting occurs for reads that overlap 100% with more than one genomic feature. 
- removal of PCR duplicates assumes that UMI sequences (2x4bp) ligated to the RNA are random. We know that they are not and this means that in libraries with low complexity reads which are not true biological duplicates are lost. We err on the side of safety by valuing specificity over sensitivity.


# A note on dependencies

- works with `bpipe = 0.9.9.3` but `0.9.9.8` introduces a behavior which breaks the forwarding on uniquely mapped reads.
- assumes `pysam >= 0.9`, may or may not work with earlier versions. 

