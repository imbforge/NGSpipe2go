# smallRNA-seq pipelines

These are pipelines developed over time for the analysis of different types of small RNAs in different experimental settings.


## piRNA profiling

This pipeline has been developed with the goal of profiling **piRNAs**, with an emphasis on Zebrafish, but it should be general enough to be appliued to other species. Some steps are specific for piRNA analyses but could be easily removed to convert this to a more general smallRNA-pipeline.

[Instructions](https://github.com/imbforge/NGSpipe2go/blob/master/README.md#preparations-to-run) on how to set-up the pipeline apply.

### Dependencies

Our systems used the module environments to manage software versions, and SLURM for queuing. Up-to-date dependencies are listed here [file](https://github.com/imbforge/NGSpipe2go/blob/master/modules/SmallRNAseq/tool.versions.groovy), but should not differ much from the following:

- Java v1.8.0_101
- R v3.5.1_debian9
- python v2.7.5
- fastqc v0.11.5
- bowtie v0.12.8
- samtools v1.5_debian9
- bedtools v2.27.1_debian9
- picard v2.7.0
- kentutils v302
- deeptools v2.4.2
- cutadapt v1.9
- fastx v0.0.14
- bowtie2 v2.3.2
- pingpongpro v1.0
- seqtk v1.2
- mirdeep2 v2.0.0.8
- repenrich v1.2
- subread v1.6.2
- htseq v0.9.0
- fastqscreen v0.12.2
- python modules:
    + pybedtools v0.7.10
    + pysam v0.13.0
- R libraries:
    + argparser
    + biomaRt
    + data.table
    + dplyr
    + edgeR
    + GeneOverlap
    + GenomicRanges
    + ggbio
    + ggplot2
    + ggthemes
    + grid
    + lemon
    + makeitprettier
    + RColorBrewer
    + scales
    + VennDiagram
- [PingPongPro](https://sourceforge.net/projects/pingpongpro/)


#### Files required
- raw reads (fastq) (though it can be changed to accept bam files)
- Genomic features (see piRNA quantification)


### What it will it do

The pipeline has several steps, some general of data preparation and some more specific for piRNA analysis. Bellow is an explanation for the most relevant, and the rationale for some of the options made. The goal is to profile small RNA-Seq libraries (RIP or input) insights into changes to piRNA pathways:
- ping-pong signal strength (z-score)
- piRNA nucleotide bias composition
- read length profiling
- Quantification of the number of small RNAs mapping to particular genomic features

The insights are limited and in most cases more in-depth analysis will be needed, but the resulting alignments and results should facilitate this. 


### read-preprocessing
We follow a very straight forward QC -> adapter removal -> low quality sequence removal strategy. Next PCR duplicates are taken care of (see bellow). There is a QC step, `FastQ screen` that looks for the presence of various contaminants in the library.


#### Filter Duplicates
As piRNAs are derived from repetitive sequences, it is expected that a large number of reads will have the same sequence. To separate what are unique RNA molecules from PCR duplicates, during the library preparation Unique Molecule Identifiers (UMIs, aka barcodes) are added to the sample before PCR amplification. These UMIs are then used to collapse the duplicated reads (PCR duplicates). Barcodes are then removed prior to mapping.

A second QC step is done at this stage (`data/qc/fastqc`), and the number of reads filtered at each stage is quantified and plotted (`results/processed_reads/figure`). 


### Mapping

Mapping to the genome is performed with bowtie with the following (relevant settings):
- **-v 2**, 2 mismatches allowed. This is allows sensitive mapping of reads in regions that might have SNPs.
- **-M 1**, if a read has more than **1** reportable alignments, one is reported at random.
- **--tryhard --best --strata --chunkmbs 256**, bowtie *best* mode.

There is also an extra step that filter the mapped files for reads which are uniquely mapped (Q 255) and keeps those separately (`results/mapped/unique`). 

The choice to retain the multimapping reads was done in the interest of flexibility - a stage has been included to separate uniquely mapping reads. The rationale is that some analysis require sensitivity, that is as many reads as possibly, for instance to quantify the classes of transposons. In this case where **exact position** of the read originated from is not very important, but rather which genomic feature is the piRNA targeting. Specificity, knowing **exactly** the genomic position of the element that piRNA is targeting is important to determine some features of phasing (Mohn *et al* 2015), for instance in the stage nucleotide signature.

For visual inspection of alignments there `bigwig` files are also created, and so are bam/bw files with separated alignments for the + and - strand.

**Ouputs:**
- `mapped/unique`
- `mapped/multimapped`
- `tracks`


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

- `PingPongPro`, which looks for positions with strong ping-pong signature (read the [manual](https://sourceforge.net/projects/pingpongpro/) for more details). By default outputs ping-pong strength per TE.
- custom script located in `NGSpipe2go/tools/piRNA/ping-pong_signature.py` that calculates the 10nt overlap bias using all input reads - by default filtered for those that map to TEs.

Whilst `PingPongPro` gives more detailed information, for instance genomic locations that could have an interesting signal, or a repeat region that might be particularly targeted, `ping-pong_signature.py` provides a quick global overview of the library to immediately identify conditions where the ping-pong cycle might be affected.

**Ouputs:**
In `results/ping-pong/`, each library will contains a plot (`figures`) (with z-score) created by `ping-pong_signature.py`, and the list of putative ping-pong signatures to the file `ping-pong_signatures.tsv` from `PingPongPro`.


#### Quantification of piRNAs

The goal of this step, `CountReads`, is to classify reads (small RNAs) according to their genomic location. That is, we try to quantify how many reads are miRNA, rRNA, repeats, etc. This way we build a profile of our small RNA library. Currently the default is to use only the TE elements annotation (DNA and RNA repeats) and plot the small RNA length profile and strand bias, but could be easily extended.

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

This format allows easy aggregation after intersecting reads with features. If the format of your `features.bed` does not conform the pipeline will fail. 


**Ouputs:**
Plots summarizing:
    - read (piRNA) length distribution
    - tables with read and matching genomic feature for further processing. 

**Limitations:**
- reasonably memory intensive because summarization is done in R. Even though `data.table` is used, it might still run into memory problems for large libraries.
- some double counting of reads but tests suggest that very few reads are double counted (if exons from protein coding genes are merged). The alternative would a stratified approach in which a read that matches more than one feature, say snoRNA and LINE, would be attributed arbitrarily to one.

### TE targeting

To determine how many piRNAs are targeting a particular TE, we make use of `RepEnrich` to count the number of reads mapping to a set of repeat elements. [RepEnrich](https://github.com/nskvir/RepEnrich) takes a two-step approach when mapping to repetitive elements:

> In RepEnrich, reads are initially aligned to the unmasked genome and divided into uniquely mapping and multi-mapping reads. Uniquely mapping reads are tested for overlap with repetitive elements, while multi-mapping reads are separately aligned to repetitive element assemblies representing individual repetitive element subfamilies (Figure 1). Repetitive element assemblies are represented by all genomic instances (assembled from the RepeatMasker annotation) of an individual repetitive element subfamily, including flanking genomic sequences, concatenated with spacer sequences to avoid spurious mapping of reads spanning multiple instances. The repetitive element assemblies are an extension of the strategy used by Day et al.[13], which however only used reads that could be unambiguously assigned to an individual subfamily.
> source: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-583


**Outputs**

In `results/repEnrich` one folder per library with table of counts for each TE family, class, or element. These can be used for differential analysis with `edgeR` or `DESeq2`.

**Limitations**

Incredibly slow. For some libraries it takes more than 24h with 4 cores. 


### Other outputs
For convenience (visualization and further custom analysis) there are several bigwigs and alignment files (bam) generated, including bigwig files for all alignments and unique (`mapped`, `/tracks/`), also separated by strand. There is also a `mappedReads.txt` file in the `mapped` folder that can be used to normalized counts to library size.  


### Summary of results

There are plans to have a report of all results being generated at the end. While this does not happen, all relevant plots are merged in a single file for a quick overview: `results/all_plots.pdf`.


## RNAi in _C. elegans_

The variety of endogenous RNAi pathways in _C. elegans_ is quite [astounding](http://wormbook.org/chapters/www_endoRNAipathwys/endoRNAipathwys.html), with a variaty of small RNA classes taking part (21U-RNA, 22G-RNA, 26G-RNA)and therefore deserves a specialized pipeline to make it justice.  

Apart from QC + read pre-processing (see piRNA pipeline above), this pipeline will:
    - summarize the library small RNa composition
    - split captured small RNAs in different classes
    - quantify these varied small RNAs


### Mapping

Mapping to the genome is performed with bowtie with the following (relevant settings):
- **-v 0**, no mismatches allowed. This is quite strict but we still get mapping rates > 80% on a regular basis. Could be changed.
- **-M 1**, if a read has more than **1** reportable alignments, one is reported at random.
- **--tryhard --best --strata --chunkmbs 256**, bowtie *best* mode.

Results are placed in `mapped/`.


### small RNA length and nucleotide composition

To have quick overview of the type of small RNAs that composed the libraries, a plot summarizing the read lengths and 5' nucleotide composition for each length. With this plot it becomes immediately apparent if a particular condition is affecting the numbers of, say 21U-RNAs. Because 26G-RNAs are usually less abundant than the other classes, a separate plot with a zoomed view of 25-28 nucleotides is also produced in a separate plot. The tables of with the counts used to produce the plots is also kept and can be used to generate more custom plots.

**outputs**

In `small_RNA_classes/sequence_bias`:

- Plots with length and 5' nucleotide composition
- tables of summaries


### How do we define small RNA classes?

Even though one could define _C. elegans_ siRNA classes strictly using the length **and** the identity of the first nucleotide **and** the mapping (sense or antisense) to particular group of target genes, _in vivo_ the size and/or nucleotide composition of each class is a little less well defined. To strike a balance between choosing the reads that are very likely to be siRNAs whilst trying to avoid excluding too many reads, our current strategy to classify reads as siRNAs is:

- 21U-RNAs (piRNAs), reads 21 nt long, and map **sense** to annotated piRNAs;
- 22G-RNAs, reads 20-23 nucleotides long, and map **antisense** to annotated genes, transposable elements or pseudogenes;
- 26G-RNAs, reads 26 nt long, and map antisense to annotated genes/pseudogenes.
- miRNA, currently not implemented.

The implementation of these rules can be found in this [stage](https://github.com/imbforge/NGSpipe2go/blob/master/modules/SmallRNAseq/filter_smRNA_classes.module.groovy). The script `NGSpipe2go/tools/smallRNA/filterSmallRNAclasses.py` can also be used to implement a stricter definition of siRNAs using also the fist base composition. 

The annotations for the gene classes can be derived from the ensemble GTFs and feed to the pipeline in the `essential.vars.groovy`. See example [here](https://github.com/imbforge/NGSpipe2go/blob/master/modules/SmallRNAseq/essential.vars.groovy#L15-18).


### Filter and quantification of small RNA classes

Reads that belong a particular siRNA class are filtered into separate bam files, which can then be used for a number of analysis. One such analysis, is the qualification of the number of reads that belong to each class. 

**output**

- bam files with the alignments that belong to each siRNA class
- plot summarizing the abundance of each class

The are located in `results/small_RNA_classes`.


### small RNA targeting

The above will tell us the global levels of siRNAs in the samples. Sometimes it is also also useful to determine which loci are affected in a condition, in the case of 21U-RNAs, or the levels of targeting in a particular gene. For that HTseq-count is used to count the number of reads per genomic feature. A few notes on this:

- HTseq-count is because it assigns reads that overlap multiple loci to that with the largest overlap (`intersection-nonempty`). This is useful because a large number of 21U-RNA loci have overlapping genomic positions. 
- by default counting uses the library preparation stranding information (`essential.vars.groovy`, variable `ESSENTIAL_STRANDED`). This means that depending on the type of small RNA of interest, 21U-RNA and miRNA should be counted sense, 22G-RNAs and 26G-RNAs should be counted antisense to features, one might have re-run this step switching the `ESSENTIAL_STRANDED` variable.

The tables of counts can be used for custom differential analysis `DESeq2` or `edgeR` analysis.   


## 21U sensor analysis

This is for a very specific type of analysis. The Ketting lab has a system in which strains containing a sensor sequence for the production of secondary 22G-RNAs, which are then crossed with mutant strains to detect 21U RNA‐induced silencing. For reference see [Fig. 3B](http://emboj.embopress.org/content/embojnl/31/16/3422/F3.large.jpg) of [Luteijn et. al.](http://emboj.embopress.org/content/31/16/3422).

In practice, in this pipeline is a follow up of the _C. elegans_ small RNA pipeline but reads are aligned to a custom reference genome that includes the sensor sequence as an extra contig. Small RNAs classes are then filtered, and the 22G-RNA coverage for the sensor is plotted. 

**outputs**

- bam files with small RNA classes
- plot of sensor coverage for each library 

## microRNA analysis

QC and read processing as described aboved (piRNA pipeline). Uses mirDeep2 to quantify miRNAs. Please use the pipeline from the [BCF](https://github.com/imbforge/NGSpipe2go/blob/master/pipelines/smallRNAseq_BCF/) instead because it is under active development. 

### Prerequisites
#### Programs required
- mirDeep2

#### Files required

- fastq.gz
- hairpin.fa provides the precursors, and
- mature.fa the mir sequences
- 
Available in iGenomes download.
