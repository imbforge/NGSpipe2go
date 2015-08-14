//rule for task bam2bw from catalog NGS, version 1
//desc: Create scaled bigwig tracks from a bam file
bam2bw = {
	doc title: "bam2bw",
		desc:  "Convert BAM file to bigWig",
		constraints: "none.",
		author: "Sergi Sayols"

	output.dir=TRACKS
	
	transform(".bam") to ("_scaled.bw") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_BEDTOOLS}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			TOTAL_MAPPED=\$( ${TOOL_SAMTOOLS} flagstat $input | head -n1| cut -f1 -d" ") &&
			SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
			genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} -g $BAM2BW_CHRSIZES > ${output.prefix}.bedgraph &&
			${TOOL_DEPENDENCIES}/ucsc/default/bedGraphToBigWig ${output.prefix}.bedgraph $BAM2BW_CHRSIZES $output
		""","bam2bw"
	}
}
