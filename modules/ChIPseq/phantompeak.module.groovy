//rule for task phantompeak from catalog ChIPseq, version 1
//desc: Phantompeak
phantompeak = {
	doc title: "Phantompeak QC  plot",
		desc:  "Phantompeak",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = QC + "/phantompeak"

	def PHANTOMPEAK_FLAGS = PHANTOMPEAK_MINSHIFT + " " + // left 'x' coordinate in plot
                            PHANTOMPEAK_MAXSHIFT + " " + // right 'x' coordinate in plot
                            PHANTOMPEAK_BINSIZE  + " " + // stepsize for cc calculation
                            PHANTOMPEAK_READLEN  + " " + // read length
                            PHANTOMPEAK_THREADS  + " " + // cores to use
                            PHANTOMPEAK_EXTRA

	transform(".bam") to("_phantompeak.png") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			source ${TOOL_R}/env.sh &&
			
			echo 'VERSION INFO'  1>&2 ; 
			echo \$(${TOOL_R}/bin/Rscript --version 2>&1 | cut -d' ' -f5) 1>&2 ;
			echo '/VERSION INFO' 1>&2 &&
			
			${TOOL_R}/bin/Rscript ${TOOL_ENCODEqc}/phantompeak.R $input \$(basename $input.prefix) $PHANTOMPEAK_FLAGS &&
			mv *_phantompeak.* $output.dir
		""","phantompeak"
	}

	forward input
}

