//vars for task Bowtie from catalog smallRNAseq_BCF, version 0.1

BOWTIE_PATH=TOOL_BOWTIE // bowtie bin to use
BOWTIE_THREADS=8	// threads to use
BOWTIE_REF=ESSENTIAL_BOWTIE_REF // prefix of the bowtie reference genome
BOWTIE_MM=1             // number of mismatches allowed
BOWTIE_MULTIREPORT=1    // if a read has more than <int> reportable alignments, one is reported at random.
BOWTIE_BEST="--tryhard --best --strata --chunkmbs 256"	// bowtie best mode
BOWTIE_QUALS="--phred33-quals"	// phred33-quals. Use --phred64-quals for old sequencing runs



