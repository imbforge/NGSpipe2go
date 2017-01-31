//rule for task extend from catalog NGS, version 1
//desc: Extend read length to the average fragment size
extend = {
	doc title: "extend",
		desc:  "Extend read length to the average fragment size",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir=MAPPED
	
	def SAMTOOLS_SORT_FLAGS = "-O bam " + EXTEND_SAMTOOLS_THREADS

	transform(".bam") to ("_ext.bam") {
		exec """
            module load bedtools/${BEDTOOLS_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&

			if [ ! -d $TMP ]; then
				mkdir -p $TMP;
			fi &&
			
			CHRSIZES=${TMP}/\$(basename ${input.prefix}).extend.chrsizes  &&
			samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			bedtools bamtobed -split -i $input |
			bedtools slop -g \$CHRSIZES -l 0 -r $EXTEND_FRAGLEN -s |
			bedtools bedtobam -ubam -g \$CHRSIZES |
			samtools sort $SAMTOOLS_SORT_FLAGS -T $TMP/\$(basename $output.prefix) - > $output &&
			samtools index $output
		""","extend"
	}
}

