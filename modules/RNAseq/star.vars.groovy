// local settings
STAR_MAXRAM = "--limitGenomeGenerateRAM 31000000000"// around 30Gb for mammals
STAR_BUFSIZE = "--limitIObufferSize 150000000"	    // buffer size
STAR_MM = "--outFilterMismatchNmax 2"               // number of mismatches allowed
STAR_MULTIMAP = "--outFilterMultimapNmax 10"        // max multimap positions per read
STAR_MININTRO = "--alignIntronMin 21"               // minimum intron size
STAR_UNMAPPED_BAM = "--outSAMunmapped Within "      // report unmapped reads to bam file? (choices: Within, None)
STAR_UNMAPPED_OUT = "--outReadsUnmapped None"       // report unmapped reads to separate fasta/fastq files, besides --outSAMunmapped? (choices: Fastx, None)


// settings imported from essential vars
STAR_THREADS = "--runThreadN " + Integer.toString(ESSENTIAL_THREADS)
STAR_REF = "--genomeDir " + ESSENTIAL_STAR_REF
STAR_OVERHANG = "--sjdbOverhang " + Integer.toString(ESSENTIAL_READLENGTH - 1)
STAR_GTF = "--sjdbGTFfile " + ESSENTIAL_GENESGTF    // gene model
STAR_SAMTOOLS_THREADS="-@ " + Integer.toString(ESSENTIAL_THREADS)

STAR_FILTER_SEC="NO"		// filter out secondary alignments from the bam file?
STAR_EXTRA=""                // extra parms to sent to the tool

// MMR settings imported from essential vars
MMR_THREADS = "--threads " + Integer.toString(ESSENTIAL_THREADS)
MMR_ANNOTATION = "--annotation " + ESSENTIAL_GENESGTF
MMR_WINDOWSIZE = "-w 60" // window size
MMR_EXTRA = ""

