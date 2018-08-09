// local settings
STAR_MM = "--outFilterMismatchNoverLmax 0.04"		// max. percentage of mismatches in the alignment (4% used by ENCODE)
STAR_UNMAPPED_BAM = "--outSAMunmapped Within"		// report unmapped reads to bam file? (choices: Within, None)

// settings imported from essential vars
STAR_THREADS = "--runThreadN " + Integer.toString(ESSENTIAL_THREADS)
STAR_REF = "--genomeDir " + ESSENTIAL_STAR_REF
STAR_OVERHANG = "--sjdbOverhang " + Integer.toString(ESSENTIAL_READLENGTH - 1)
STAR_GTF = "--sjdbGTFfile " + ESSENTIAL_GENESGTF		// gene model GTF file

// additional settings
STAR_EXTRA=""		// extra parms to sent to the tool
STAR_FILTER_SEC="YES"		// remove secondary alignments in order to keep just 1 alignment of the multi-mapping reads ?
STAR_SAMTOOLS_THREADS= "-@ " + Integer.toString(ESSENTIAL_THREADS)
