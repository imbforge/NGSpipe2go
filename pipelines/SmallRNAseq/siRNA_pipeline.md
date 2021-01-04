# siRNA and siRNA sensor pipelines

This pipeline is specifically to analysis siRNAs in _C elegans_ (worm), but could be applied to other _Caenorhabditis_ species. This is becasue siRNA mechanisms of this taxa is so complex and different from other organisms. 

The `siRNA_sensor` pipeline is basically the same but adds stages to anlysis siRNA activity in a 21U sensor. See [Placentino et al.](https://www.embopress.org/doi/full/10.15252/embj.2020105280)  

Brief outlook of the key outputs:

- read length profile which is a good proxy to establish the main classes of siRNAs in the libraries
- generate bam files containing each of the siRNA classes (21U, 22G, 26G, miRNA)
- read count siRNas mapping to genomic features

Just like siRNA pathways in _C. elegans_ are complex, so are the biological questions and follow up analysis which can be done on sequencing data. Therefore, this pipeline is geared towards generating intermediate files to be re-used for custom analysis.

# Inputs

- raw reads
- genome sequence (fasta)
- genome annotation (gtf)
- subset GTFs needed to filter siRNA classes.
- rRNA sequences

For examples of how to set-up the files see the [test data folder](testdata/siRNA/README.md)


# Main steps and parameters

- QC
    + adapter trimming
    + removal of reads with low-quality bases
    + removal of PCR duplicates
    + FastQC, performed before and after mapping
    + diagnostics plots
- summarize read lengths
- mapping with bowtie, allowing for one randomly selected alignment for multi-mapping reads, and no mismatches
- Filter non-structural reads (non rRNA/tRNA/snoRNA/snRNA). Only these will be use for further analysis.
- separate 21U, 22G, 26G and miRNA reads. See below for criteria used. Can be changed in the file `modules/SmallRNAseq/filter_smRNA_classes.vars.groovy`
- count reads with HTseq-count. This was used because will allows to disambiguate reads that partially overlap more than one genomic feature. Very important for 21Us whose coding units overlap in the genome.
- generate bigwigs 

The **sensor pipeline** contains in addition:

- mapping is done to a genome file contained the sensor sequence in addition to the canonical worm chromosome sequences
- Filter reads mapping to the sensor
- plot 22G coverage of sensor


# siRNA definition

- `21U` are 21 nt long sequences mapping sense to annotated piRNA loci and starting with a T;
- `22G` are 20-23 nt, starting with a G, and map antisense to protein-coding/pseudogenes/lincRNA/transposons;
- `26G`, are those which are 26 nt, starting with a G, and map antisense to annotated protein-coding/pseudogenes/lincRNA; 
- `miRNAs` are 20-24 nt mapping sense to annotated miRNA loci

# Limitations and future improvements

- It could have more downstream analysis stages which are fairly common to most projects.
- removal of PCR duplicates assumes that UMI sequences (2x4bp) ligated to the RNA are random. We know that they are not and this means that in libraries with low complexity reads which are not true biological duplicates are lost. We err on the side of safety by valuing specificity over sensitivity.
