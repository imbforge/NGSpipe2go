//rule for task Bam2bw from catalog NGS, version 1
//desc: Create scaled bigwig tracks from a bam file
BamCoverageStrands = {
	doc title: "BamCoverage",
		desc:  "Convert BAM file to bigWig. The --filterRNAstrand option assumes the sequencing library generated from ILLUMINA dUTP/NSR/NNSR methods",
		constraints: "Requires deepTools >= 2.2",
		author: "Antonio Domingues"

	output.dir=TRACKS

	transform(".bam") to (".fwd.scaled.bw", ".rev.scaled.bw") {
     	exec """
        module load bedtools/${DEEPTOOLS_VERSION} &&

        bamCoverage --bam $input -o $output1 \
            --normalizeTo1x $GENOME_SIZE \
            --binSize 10 \
            --numberOfProcessors $ESSENTIAL_CORES \
            --filterRNAstrand forward &&

        bamCoverage --bam $input -o $output1 \
            --normalizeTo1x $GENOME_SIZE \
            --binSize 10 \
            --numberOfProcessors $ESSENTIAL_CORES \
            --filterRNAstrand reverse 
            ""","BamCoverageStrands"
	}
}

