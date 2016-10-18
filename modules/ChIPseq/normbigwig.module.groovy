normbigwig = {
	doc title: "normbigwig",
		desc:  "normbigwig wrapper",
		constraints: "Only performs treatment control substracted normlised bigwig tracks",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Nastasja Kreim"

	output.dir = NORMBIGWIG_OUTDIR 

	transform(".bam") to("_normbigwig.done") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			source ${TOOL_DEEPTOOLS}/env.sh &&
			source ${TOOL_R}/env.sh &&
			source ${TOOL_UCSC}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi;
			if [ ! -d ${TMP} ]; then
				mkdir -p ${TMP};
			fi &&
			if [ ! -e $NORMBIGWIG_TARGETS ]; then
				echo "Targets file $NORMBIGWIG_TARGETS doesn't exist" >> $output &&
				exit 0;
			fi;
			touch $output;
			BAM=\$(basename $input) &&
			grep \$BAM $NORMBIGWIG_TARGETS | while read -r TARGET; do
				IP=\$(       echo \$TARGET | cut -f1 -d" ") &&
				IPname=\$(   echo \$TARGET | cut -f2 -d" ") &&
				INPUT=\$(    echo \$TARGET | cut -f3 -d" ") &&
				INPUTname=\$(echo \$TARGET | cut -f4 -d" ");
				if [ "\$BAM" != "\$INPUT" ]; then
			echo "\${IPname} vs \${INPUTname}" >> $output ;
			CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
			${TOOL_SAMTOOLS} idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			bamCompare -b1 $input -b2 $NORMBIGWIG_MAPPED/$INPUT --numberOfProcessors $NORMBIGWIG_THREADS $NORMBIGWIG_OTHER --outFileName $TMPDIR/\${BAM%.bam}_\${INPUTname}_norm.bedgraph &&
			sort -k1,1 -k2,2n  $TMPDIR/\${BAM%.bam}_\${INPUTname}_norm.bedgraph >  $TMPDIR/\${BAM%.bam}_\${INPUTname}_norm.bedgraph.sorted && 
			bedGraphToBigWig \${TMPDIR}/\${BAM%.bam}_\${INPUTname}_norm.bedgraph.sorted $CHRSIZES  $output.dir/\${BAM%.bam}_\${INPUTname}_norm.bw;
			if [ \$? -ne 0 ]; then rm $output; fi;
				fi;
			done
		""","normbigwig"
	}
	forward input
}

