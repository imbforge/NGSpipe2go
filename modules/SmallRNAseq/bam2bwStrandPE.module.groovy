Bam2bwStrandPE = {
	doc title: "Bam2bw",
		desc:  "Convert BAM file to bigWig for each strand. Normalization factor is the total number of reads. First separates reads. Source and inspiration: https://www.biostars.org/p/92935/",
		constraints: "none.",
		author: "Antonio Domingues"

	def EXP = input.split("/")[-1].replaceAll(".bam", "")
	output.dir=BAMCOVSTRANDSPE_OUTDIR

	transform(".bam") to (".fwd_pe.scaled.bw", ".rev_pe.scaled.bw")  {
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

			samtools view -b -f 128 -F 16 $input > ${TMP}/${EXP}.fwd1.bam &&
			samtools index ${TMP}/${EXP}.fwd1.bam &&

			samtools view -b -f 80 $input > ${TMP}/${EXP}.fwd2.bam &&
			samtools index ${TMP}/${EXP}.fwd2.bam &&

			samtools merge -f ${TMP}/${EXP}.fwd.bam ${TMP}/${EXP}.fwd1.bam ${TMP}/${EXP}.fwd2.bam &&
			samtools index ${TMP}/${EXP}.fwd.bam &&

			samtools view -b -f 144 $input > ${TMP}/${EXP}.rev1.bam &&
			samtools index ${TMP}/${EXP}.rev1.bam &&

			samtools view -b -f 64 -F 16 $input > ${TMP}/${EXP}.rev2.bam &&
			samtools index ${TMP}/${EXP}.rev2.bam &&
	
			samtools merge -f ${TMP}/${EXP}.rev.bam ${TMP}/${EXP}.rev1.bam ${TMP}/${EXP}.rev2.bam &&
			samtools index ${TMP}/${EXP}.rev.bam &&

			genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${TMP}/${EXP}.fwd.bam -g \${CHRSIZES} > ${output1.prefix}.bedgraph &&
			bedGraphToBigWig ${output1.prefix}.bedgraph \${CHRSIZES} $output1 &&
			
			genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${TMP}/${EXP}.rev.bam -g \${CHRSIZES} > ${output2.prefix}.bedgraph &&
			bedGraphToBigWig ${output2.prefix}.bedgraph \${CHRSIZES} $output2 &&

			rm ${CHRSIZES} ${output1.prefix}.bedgraph ${output2.prefix}.bedgraph ${TMP}/${EXP}.*.bam 

		""","Bam2bwStrandPE"
	}
}

