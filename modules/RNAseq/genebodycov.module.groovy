//rule for task geneBodyCov from catalog RNAseq, version 1
//desc: Calculate the RNA-seq coverage over gene body
geneBodyCov = {
	doc title: "geneBodyCoverage",
		desc:  """Calculate the RNA-seq coverage over gene body. 
			Useful to check the 5' or 3' coverage bias""",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

    output.dir = GENEBODYCOV_OUTDIR
	
	// println("DEBUG: " + output.dir)

    // run the chunk
    //from(".bam") produce(input.prefix + ".geneBodyCoverage.curves.png",
    //                     input.prefix + ".geneBodyCoverage.r",
    //                     input.prefix + ".geneBodyCoverage.txt") {
	transform(".bam") to (".geneBodyCoverage.curves.png", ".geneBodyCoverage.r", ".geneBodyCoverage.txt") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES && 
			source ${TOOL_RSeQC}/env.sh && 
			if [ -n "\$LSB_JOBID" ]; then 
				export TMPDIR=/jobdir/\${LSB_JOBID}; 
			fi &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(python ${TOOL_RSeQC}/bin/geneBody_coverage.py --version) 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			
			python ${TOOL_RSeQC}/bin/geneBody_coverage.py -i $input -f png -r $GENEBODYCOV_BED -o ${output3.prefix.prefix}
		""","geneBodyCov"
	}
	// python ${TOOL_RSeQC}/bin/geneBody_coverage.py -i $input -f png -r $GENEBODYCOV_BED -o ${output.prefix} && mv ${input.prefix}.geneBodyCoverage.* $output.dir

	forward input

}
