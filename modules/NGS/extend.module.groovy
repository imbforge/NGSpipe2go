//rule for task extend from catalog NGS, version 1
//desc: Extend read length to the average fragment size
extend = {
	doc title: "extend",
		desc:  "Extend read length to the average fragment size",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir=MAPPED

	transform(".bam") to ("_ext.bam") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_BEDTOOLS}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&
			
			if [ ! -d ${TMP} ]; then
				mkdir -p ${TMP};
			fi &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(bedtools --version |  cut -d' ' -f2) 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			
			CHRSIZES=${TMP}/\$(basename ${input.prefix}).extend.chrsizes  &&
			${TOOL_SAMTOOLS} idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			bedtools bamtobed -split -i $input |
			bedtools slop -g \$CHRSIZES -l 0 -r $EXTEND_FRAGLEN -s |
			bedtools bedtobam -ubam -g \$CHRSIZES |
			${TOOL_SAMTOOLS} sort -@ $EXTEND_SAMTOOLS_THREADS -T $TMP/\$(basename $output.prefix) - > $output &&
			${TOOL_SAMTOOLS} index $output
		""","extend"
	}
}

