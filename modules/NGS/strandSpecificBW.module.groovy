strandBigWig = {
	doc title: "strandBigWig",
		desc:  "strandBigWig wrapper",
		constraints: "strandspecific bigwig for rnaseq",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Nastasja Kreim"

	output.dir = STRANDSPECIFICBIGWIG_OUTDIR
	//this might be confusing regarding the reverse and forward
	//setting but according to the deeptools manual it has to be like
	//that. 
    	if(ESSENTIAL_STRANDED == "yes") {
		FORWARD="reverse"
		REVERSE="forward"
	} else if(ESSENTIAL_STRANDED == "reverse") {
		FORWARD="forward"
		REVERSE="reverse"
	}

	transform(".bam") to(".fwd.bw", ".rev.bw") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			module load deepTools/$DEEPTOOLS_VERSION &&
			module load samtools/$SAMTOOLS_VERSION &&
			module load kentUtils/$KENTUTILS_VERSION &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi;
			base=\$(basename $input.prefix) &&
			echo \$base &&
			CHRSIZES=${TMPDIR}/\${base}.bam2bw.chrsizes &&
			samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
			bamCoverage --numberOfProcessors $STRANDSPECIFICBIGWIG_CORES --filterRNAstrand $FORWARD $STRANDSPECIFICBIGWIG_OTHER -b  $input -o \${TMPDIR}/\${base}.fwd.bedgraph &&
			bamCoverage --numberOfProcessors $STRANDSPECIFICBIGWIG_CORES --filterRNAstrand $REVERSE $STRANDSPECIFICBIGWIG_OTHER -b $input -o \${TMPDIR}/\${base}.rev.bedgraph&&
			
			awk 'BEGIN {OFS="\t"; FS="\t"}{print \$1,\$2, \$3,"-"\$4}' \${TMPDIR}/\${base}.rev.bedgraph > \${TMPDIR}/\${base}.rev.bedgraph_neg &&
			mv $TMPDIR/\${base}.rev.bedgraph_neg $TMPDIR/\${base}.rev.bedgraph &&
			bedGraphToBigWig \${TMPDIR}/\${base}.fwd.bedgraph $CHRSIZES $output1 &&
			bedGraphToBigWig \${TMPDIR}/\${base}.rev.bedgraph $CHRSIZES $output2
			
		""","strandBigWig"
	}
	forward input
}


