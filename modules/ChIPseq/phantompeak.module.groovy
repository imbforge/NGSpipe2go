//rule for task phantompeak from catalog ChIPseq, version 1
//desc: Phantompeak
phantompeak = {
	doc title: "Phantompeak QC  plot",
		desc:  "Phantompeak",
		constraints: "",
		author: "Sergi Sayols"

	output.dir = QC + "/phantompeak"

	def PHANTOMPEAK_FLAGS = Integer.toString(PHANTOMPEAK_MINSHIFT) + " " + // left 'x' coordinate in plot
	                        Integer.toString(PHANTOMPEAK_MAXSHIFT) + " " + // right 'x' coordinate in plot
							Integer.toString(PHANTOMPEAK_BINSIZE)  + " " + // stepsize for cc calculation
							Integer.toString(PHANTOMPEAK_READLEN)  + " " + // read length
							Integer.toString(PHANTOMPEAK_THREADS)          // cores to use

	transform(".bam") to("_phantompeak.png") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi;
			
			echo 'VERSION INFO'  1>&2 &&
			${TOOL_R}/bin/Rscript --version 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&
			
			${TOOL_R}/bin/Rscript ${TOOL_ENCODEqc}/phantompeak.R $input \$(basename $input.prefix) $PHANTOMPEAK_FLAGS &&
			mv *_phantompeak.* $output.dir
		""","phantompeak"
	}

	forward input
}

