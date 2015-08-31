//rule for task extend from catalog NGS, version 1
//desc: Extend read length to the average fragment size
extend = {
	doc title: "extend",
		desc:  "Extend read length to the average fragment size",
		constraints: "",
		author: "Sergi Sayols"

	output.dir=MAPPED

	transform(".bam") to ("_ext.bam") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_BEDTOOLS}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&

			CHRSIZES=${TMP}/\$(basename ${input.prefix}).extend.chrsizes  &&
			${TOOL_SAMTOOLS} idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			bedtools bamtobed -split -i $input |
			bedtools slop -g \$CHRSIZES -l 0 -r $EXTEND_FRAGLEN -s |
			bedtools bedtobam -ubam -g \$CHRSIZES |
			${TOOL_SAMTOOLS} sort - $output.prefix &&
			${TOOL_SAMTOOLS} index $output
		""","extend"
	}
}

