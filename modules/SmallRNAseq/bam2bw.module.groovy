//rule for task Bam2bw from catalog NGS, version 1
//desc: Create scaled bigwig tracks from a bam file
Bam2bw = {
	doc title: "Bam2bw",
		desc:  "Convert BAM file to bigWig",
		constraints: "none.",
		author: "Sergi Sayols"

	output.dir=TRACKS

	transform(".bam") to (".scaled.bw") {
		exec """
			if [ ! -d ${TMP} ]; then
				mkdir -p ${TMP};
			fi &&

			module load bedtools/${BEDTOOLS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load kentUtils/${KENTUTILS_VERSION} &&


			CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
			samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			TOTAL_MAPPED=\$( samtools flagstat $input | head -n1 | cut -f1 -d" ") &&
			SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
			genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} -g \${CHRSIZES} > ${output.prefix}.bedgraph &&
			bedGraphToBigWig ${output.prefix}.bedgraph \${CHRSIZES} $output &&
			rm \${CHRSIZES} ${output.prefix}.bedgraph
			
		""","Bam2bw"
	}
}

