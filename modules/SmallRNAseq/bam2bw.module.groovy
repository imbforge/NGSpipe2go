//rule for task Bam2bw from catalog NGS, version 1
//desc: Create scaled bigwig tracks from a bam file
Bam2bw = {
	doc title: "Bam2bw",
		desc:  "Convert BAM file to bigWig",
		constraints: "none.",
		author: "Sergi Sayols, Antonio Domingues"

	output.dir=TRACKS
	
	def SAMPLE_NAME = input.split("/")[-1].split("\\.", 2)[0]

	transform(".bam") to (".scaled.bw") {
		exec """
			module load bedtools/${BEDTOOLS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load kentUtils/${KENTUTILS_VERSION} &&

			if [ ! -d ${TMP} ]; then
				mkdir -p ${TMP};
			fi &&
			
			CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
			samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			
			TOTAL_MAPPED=\$( samtools flagstat $input | head -n5 | tail -n1 | cut -f1 -d" ") &&

			SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
			genomeCoverageBed -ibam -split -scale \${SCALE} -bg -i ${input} | bedSort stdin stdout > ${output.prefix}.bedgraph &&
			bedGraphToBigWig ${output.prefix}.bedgraph \${CHRSIZES} $output &&
			rm \${CHRSIZES} ${output.prefix}.bedgraph
		""","Bam2bw"
	}
}

