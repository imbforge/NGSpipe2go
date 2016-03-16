//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
BaseRecalibration = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"

	output.dir = MAPPED
    def GATK_FLAGS = "-knownSites latest_dbsnp.vcf "
    
    // check if a region limit was provided
    if (GATK_CALL_REGION!=null && GATK_CALL_REGION.length()>0) {
        GATK_FLAGS = GATK_FLAGS + " -L " + GATK_CALL_REGION
    } else {
        GATK_FLAGS = ""
    }
        
    transform (".bam") to (".recalibration.table", ".recalibrated.bam") {
        exec """
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
            echo 'VERSION INFO'  1>&2 &&
            echo \$(java -jar $TOOL_GATK/GenomeAnalysisTK.jar --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&
            
            java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -nct $GATK_THREADS -R $GATK_BWA_REF -knownSites ${GATK_KNOWN_VARIANTS} -I $input -o $output1 &&
            
            java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T PrintReads -nct $GATK_THREADS -R $GATK_BWA_REF -I $input -BQSR $output1 -o $output2 
        ""","BaseRecalibration"
    }
    
}
