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
	def GENEBODYCOV_FLAGS = GENEBODYCOV_FORMAT + " " +
                            GENEBODYCOV_BED    + " " +
                            GENEBODYCOV_EXTRA
	
    // run the chunk
	transform(".bam") to (".geneBodyCoverage.curves.png", ".geneBodyCoverage.r", ".geneBodyCoverage.txt") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES && 
			source ${TOOL_RSeQC}/env.sh && 
			if [ -n "\$LSB_JOBID" ]; then 
				export TMPDIR=/jobdir/\${LSB_JOBID}; 
			fi &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(python ${TOOL_RSeQC}/bin/geneBody_coverage.py --version | cut -d' ' -f2) 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
			
			python ${TOOL_RSeQC}/bin/geneBody_coverage.py -i $input -o ${output3.prefix.prefix} $GENEBODYCOV_FLAGS
		""","geneBodyCov"
	}
	forward input
}
