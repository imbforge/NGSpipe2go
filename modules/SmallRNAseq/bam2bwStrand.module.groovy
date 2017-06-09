Bam2bwStrand = {
	doc title: "Bam2bw",
		desc:  "Convert BAM file to bigWig for each strand. Normalization factor is the total number of reads.",
		constraints: "none.",
		author: "Antonio Domingues"

	output.dir=BAMCOVSTRANDS_OUTDIR

	transform(".bam") to (".scaled.fwd.bw", ".scaled.rev.bw")  {
		exec """
			module load bedtools/${BEDTOOLS_VERSION} &&
			module load samtools/${SAMTOOLS_VERSION} &&
			module load kentUtils/${KENTUTILS_VERSION} &&

			if [ ! -d ${TMP} ]; then
				mkdir -p ${TMP};
			fi &&
			
			CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
			samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			TOTAL_MAPPED=\$( samtools flagstat $input | head -n1| cut -f1 -d" ") &&
			SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&

			genomeCoverageBed -bg -split -scale \${SCALE} -strand "+" -ibam ${input} -g ${CHRSIZES} > ${output1.prefix}.bedgraph &&
			bedGraphToBigWig ${output1.prefix}.bedgraph ${CHRSIZES} $output1 &&
			
			genomeCoverageBed -bg -split -scale ${SCALE} -strand "-" -ibam ${input} -g ${CHRSIZES} > ${output2.prefix}.bedgraph &&
			bedGraphToBigWig ${output2.prefix}.bedgraph ${CHRSIZES} $output2 &&
			
			rm ${CHRSIZES} ${output1.prefix}.bedgraph ${output2.prefix}.bedgraph

		""","Bam2bwStrand"
	}
}

