Cutadapt = {
		doc title: "Cutadapt",
		desc:  "remove adapter from reads",
		constraints: "Only supports compressed FASTQ files",
		author: "Antonio Domingues, Anke Busch"

	output.dir = CUTADAPT_OUTDIR

	// create the log folder if it doesn't exists
	def CUTADAPT_LOGDIR = new File( CUTADAPT_OUTDIR + "/logs")
	if (!CUTADAPT_LOGDIR.exists()) {
		CUTADAPT_LOGDIR.mkdirs()
	}	

	transform(".fastq.gz") to (".cutadapt.fastq.gz",".cutadapt_discarded.fastq.gz") {

		 def SAMPLENAME = input.prefix.prefix
	
     		 exec """

            module load cutadapt/${CUTADAPT_VERSION} &&
        
			SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
			
			cutadapt $ADAPTER_SEQUENCE -O $MINIMUM_OVERLAP -m $MINIMUM_LENGTH_KEEP -M $MAXIMUM_LENGTH_KEEP -o $output1 $input 2>&1 >> ${CUTADAPT_LOGDIR}/\${SAMPLENAME_BASE}.cutadapt.log &&
			cutadapt $ADAPTER_SEQUENCE -O $MINIMUM_OVERLAP -m $MAXIMUM_LENGTH_KEEP_PLUS1 -o $output2 $input 2>&1 >> ${CUTADAPT_LOGDIR}/\${SAMPLENAME_BASE}.cutadapt_discarded.log
			
		""","Cutadapt"
	}
}
