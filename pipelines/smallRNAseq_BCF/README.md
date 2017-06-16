# smallRNA-seq pipelines

## piRNA profiling

This pipeline has been developed with the goal of profiling **piRNAs**, so keep that in mind when using it. Some steps are specific for piRNA analyses but could be easily removed to convert this to a more general smallRNA-pipeline. It is also geared to RIP-seq data. [Instructions](https://github.com/imbforge/NGSpipe2go/blob/master/README.md#preparations-to-run) on how to set-up the pipeline apply.

### Prerequisites

#### Programs required
- FastQC
- Samtools (> 1.3)
- BED tools (> 2.25.0)
- UCSC tool set
- cutadapt
- FASTX-Toolkit
- bowtie (= 0.12.8)
- python modules:
    - pybedtools
    - pysam
- R libraries:
    - ggplot2
    - RColorBrewer
    - scales
    - data.table
    - plyr
- [PingPongPro](https://sourceforge.net/projects/pingpongpro/)

#### Files required
- raw reads (fastq)
- Genomic features (see piRNA quantification)

### What it will it do
The pipeline has several steps, some general of data preparation and some more specific for piRNA analysis. Bellow is an explanation for the most relevant, and the rationale for some of the options made. The goal is to profile small RNA-Seq libraries (RIP or input) insights into changes to piRNA pathways:
- ping-pong signal strength
- piRNA nucleotide bias composition
- read length profiling
- target site bias

The insights are limited and in most cases more in-depth analysis will be needed, but the resulting alignments and results should facilitate this. It is a work in progress.

### read-preprocessing
We follow a very straight forward QC -> adapter removal -> low quality sequence removal strategy. Next PCR duplicates are taken care of (see bellow).

#### Filter Duplicates
As piRNAs are derived from repetitive sequences, it is expected that a large number of reads will have the same sequence. To separate what are unique RNA molecules from PCR duplicates, during the library preparation Unique Molecule Identifiers (UMIs, aka barcodes) are added to the sample before PCR amplification. These UMIs are then used to collapse the duplicated reads (PCR duplicates).

NOTE: at the end of this stage barcodes are still present in the read sequence and will be removed during mapping. In the future a stage might be included to remove these prior to mapping.

A second QC step is done at this stage, and the number of reads filtered at each stage is quantified and plotted (link_to_folder). 

An overview of filtered reads (PCR duplicates) can be found in `results/processed_reads/figure/`

### Mapping
Mapping to the genome is performed with bowtie with the following (relevant settings):
- **-v 0**, no mismatches allowed. Considering changes to **1** in the future.
- **-M 1**, if a read has more than **1** reportable alignments, one is reported at random.
- **--tryhard --best --strata --chunkmbs 256**, bowtie *best* mode.
- **--trim5 4 --trim3 4** to remove the barcodes (UMIs).

The choice to retain the multimapping reads was done in the interest of flexibility - a stage has been included to separate uniquely mapping reads. The rationale is that some analysis require specificity, that is as many reads as possibly, for instance to quantify the classes of transposons. In this case where exacly the read originated from is not very important, but rather which genomic feature is the piRNA targeting. Specificity, knowing **exactly** the genomic position of the element that piRNA is targeting is important to determine some features of phasing (Mohn *et al* 2015), for instance in the stage nucleotide signature.

**Ouputs:**
- `results/mapped/unique`
- `results/mapped/multimapped`

Note that uniquely mapped reads are defined by their mapping score `255`.

### piRNA analysis
#### Nucleotide signature
Determines the nucleotide bias of small RNA reads. It is expected that libraries enriched in piRNAs will have be enriched in U(ridine) at position 1 (Ziwi bound) and/or A(denine) at position 10 (Zili bound).

The nucleotide bias/signature of piRNAs is determined using the custom script `NGSpipe2go/tools/piRNA/piRNABaseTerminalBases.py`. Sense and antisense is determined using a features bed file.

**Ouput**
Can be found in `results/nucleotide_signature/` and each library in a separate folder. Results are tables with nucleotide frequency at each position, in the sense and antisense strand, 5' and 3'. It also includes +/- 20 nucleotides (default) upstream and downstream of the piRNA which Useful to see phasing - **U** bias at +1 downstream of the piRNA sequence (Mohn *et al* 2015). In the `figures` folder there are the resulting plots.

#### ping-pong signal
From PingPongPro:
>The ping-pong cycle produces short RNA molecules which are complementary at their 5’ ends over a length of exactly ten nucleotides. In piRNA-Seq data, these molecules appear as stacks of reads on opposite strands, which overlap by ten nucleotides

This means that the stronger the cycle is, stronger the signal at position 10. The pipeline calculates the ping-pong-signal using two different tools:
- `PingPongPro`, which looks for position with strong ping-pong signature (read the [manual](https://sourceforge.net/projects/pingpongpro/) for more details);
- custom script located in `NGSpipe2go/tools/piRNA/ping-pong_signature.py` that calculates the 10nt overlap bias in the full library.

Whilst `PingPongPro` gives more detailed information, for instance genomic locations that could have an interesting signal, or a repeat region that might be particularly targeted, `ping-pong_signature.py` provides a quick global overview of the library to immediately identify conditions where the ping-pong cycle might be affected.

**Ouputs:**
In `results/ping-pong/`, each library will contains a plot (`figures`) (with z-score) created by `ping-pong_signature.py`, and the list of putative ping-pong signatures to the file `ping-pong_signatures.tsv` from `PingPongPro`.

#### Quantification of piRNAs
In step reads are classified according to their genomic location. That is, we try to quantify how many reads are miRNA, rRNA, repeats, etc. This way we build a profile of our small RNA library.

Genomic features are in bed6 format, with an extra column (col7) indicating the biotype. Example:

> chr1    451140  451218  ENSDARG00000081781      .       +       miRNA
> chr1    1275347 1275428 ENSDARG00000080051      .       +       miRNA
> chr1    2689142 2689237 ENSDARG00000080082      .       -       miRNA
> chr1    2806070 2806208 ENSDARG00000082515      .       +       miRNA
> chr1    2806231 2806314 ENSDARG00000083333      .       +       miRNA
> chr1    2806383 2806472 ENSDARG00000081650      .       +       miRNA
> chr1    2806564 2806717 ENSDARG00000082292      .       +       miRNA
> chr1    2806749 2806836 ENSDARG00000081789      .       +       miRNA
> chr1    2806868 2806953 ENSDARG00000083147      .       +       miRNA
> chr1    12169738        12169820        ENSDARG00000091084      .       +       miRNA

Repeats are taken from the repeatMasked track (UCSC) and thus have repName (col4) and repFamily (col5). repClass is col7 (~ biotype):

> chr1    76692   76762   Rex-Babar       CR1-31_DR       -       LINE
> chr1    85496   85570   Rex-Babar       CR1-31_DR       +       LINE
> chr1    86642   86716   Rex-Babar       CR1-31_DR       -       LINE
> chr1    115656  115790  L2      CR1-1_DR        -       LINE
> chr1    118919  119087  RTE     EXPANDER1_DR    -       LINE
> chr1    137923  138003  Rex-Babar       CR1-40_DR       +       LINE
> chr1    141503  141590  Rex-Babar       CR1-26_DR       +       LINE
> chr1    146719  146803  I       I-2_DR  +       LINE
> chr1    151889  152160  I       I-3_DR  +       LINE
> chr1    176929  177671  Rex-Babar       CR1-26_DR       -       LINE 

This format allows easy aggregation after intersecting reads with features. If the format of your `features.bed` is not this the pipeline will fail. Any arbitrary number/type of biotypes with the exception of repeat classes that need to be present.

**Ouputs:**
Plots summarizing:
- read (piRNA) length distribution,
- relative abundance of reads mapping to biotypes;
- distribution of piRNAs per repeat class;

**Limitations:**
- reasonably memory intensive because summarization is done in R. Even though `data.table` is used, it might still run into memory problems for large libraries.
- some double counting of reads but tests suggest that very few reads are double counted (if exons from protein coding genes are merged). The alternative would a stratified approach in which a read that matches more than one feature, say snoRNA and LINE, would be attributed arbitrarily to one.


### Other outputs
For convenience (visualization and further custom analysis) there are several bigwigs and alignment files (bam) generated, including bigwig files for all alignments and unique (`results/mapped/tracks/`), also separated by strand.


### Summary of the pipeline
There are plans to have a report of all results being generated at the end. While this does not happen, all relevant plots are merged in a single file for a quick overview: `results/all_plots.pdf`.


## microRNA analysis

To be properly documented. Uses mirDeep2 to quantify miRNAs.

### Prerequisites
#### Programs required
- mirDeep2

#### Files required
- hairpin.fa provides the precursors, and
- mature.fa the mir sequences
Available in iGenomes download.
