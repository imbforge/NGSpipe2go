//rule for task bam2bw from catalog NGS, version 1
//desc: Create scaled bigwig tracks from a bam file
bam2bw = {
	doc title: "bam2bw",
		desc:  "Convert BAM file to bigWig",
		constraints: "none.",
		author: "Sergi Sayols"

	output.dir=TRACKS

	transform(".bam") to (".scaled.bw") {
		exec """
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			if [ ! -d ${TMP} ]; then
				mkdir -p ${TMP};
			fi &&

			echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_BEDTOOLS}/genomeCoverageBed -h 2>&1 | grep 'Version') 1>&2
			${TOOL_KENTUTILS}/bedGraphToBigWig    1>&2 &&
			echo '/VERSION INFO' 1>&2 &&

			CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
			${TOOL_SAMTOOLS} idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			TOTAL_MAPPED=\$( ${TOOL_SAMTOOLS} flagstat $input | head -n1 | cut -f1 -d" ") &&
			SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
			${TOOL_BEDTOOLS}/genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} -g \${CHRSIZES} > ${output.prefix}.bedgraph &&
			${TOOL_KENTUTILS}/bedGraphToBigWig ${output.prefix}.bedgraph \${CHRSIZES} $output &&
			rm \${CHRSIZES} ${output.prefix}.bedgraph
		""","bam2bw"
	}
}

