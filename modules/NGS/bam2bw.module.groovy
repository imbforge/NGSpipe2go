//rule for task bam2bw from catalog NGS, version 1
//desc: Create scaled bigwig tracks from a bam file
bam2bw = {
	doc title: "bam2bw",
		desc:  "Convert BAM file to bigWig",
		constraints: "none.",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir=TRACKS
	
	transform(".bam") to ("_scaled.bw") {
		exec """
			module load bedtools/${BEDTOOLS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load kentUtils/${KENTUTILS_VERSION} &&

			if [ -n "\$SLURM_JOBID" ]; then
                                export TMPDIR=/jobdir/\${SLURM_JOBID};
                        fi &&

			BASEOUTPUT=`basename $output` &&

			CHRSIZES=\${TMPDIR}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
                        samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
                        TOTAL_MAPPED=\$( samtools flagstat $input | head -n5 | tail -n1 | cut -f1 -d" ") &&
                        SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
                        genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} | sortBed -i -  > \${TMPDIR}/\${BASEOUTPUT%.bw}.bedgraph &&
                        bedGraphToBigWig \${TMPDIR}/\${BASEOUTPUT%.bw}.bedgraph \${CHRSIZES} \${TMPDIR}/\${BASEOUTPUT} &&
                        cp \${TMPDIR}/\${BASEOUTPUT} $output
			
		""","bam2bw"
	}
}

