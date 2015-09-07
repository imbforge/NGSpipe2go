//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
BaseRecalibration = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	author: "Oliver Drechsel"

	output.dir = MAPPED
    
    transform (".bam") to (".recalibration.table", ".recalibrated.bam") {
        
        def GATK_FLAGS = "-knownSites latest_dbsnp.vcf "
        
        
        exec """
        
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
        java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -nct $GATK_THREADS -R $ESSENTIAL_BWA_REF -knownSites ${ESSENTIAL_KNOWN_VARIANTS} -I $input -o $output1 &&
        
        java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T PrintReads -nct $GATK_THREADS -R $ESSENTIAL_BWA_REF -I $input -BQSR $output1 -o $output2 
        
        ""","BaseRecalibration"
    }
    
}