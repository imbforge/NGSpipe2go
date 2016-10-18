//rule for task infer_experiment from catalog RNAseq, version 1
inferexperiment = {
	doc title: "inferexperiment",
		desc:  """Calculate the strand-specificity. 
			of the library. """,
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Nastasja Kreim"

    output.dir = INFEREXPERIMENT_OUTDIR

    // run the chunk
    transform(".bam") to (input.prefix + "_inferexperiment.txt") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES && 
			source ${TOOL_RSeQC}/env.sh && 
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(python ${TOOL_RSeQC}/bin/infer_experiment.py --version | cut -d' ' -f2) 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			
			python ${TOOL_RSeQC}/bin/infer_experiment.py -i $input $INFEREXPERIMENT_EXTRA $INFEREXPERIMENT_BED > $output
		""","inferexperiment"
	}

	forward input
}

