//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
VariantCallHC = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"

	output.dir = RESULTS
    
    def GATK_FLAGS = ""
    
    // check if a region limit was provided
    if (ESSENTIAL_CALL_REGION!=null && ESSENTIAL_CALL_REGION.length()>0) {
        
        GATK_FLAGS = " -L " + ESSENTIAL_CALL_REGION
        
    } else {
        
        GATK_FLAGS = ""
        
    }
    
    transform (".duprm.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        exec """
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
            echo 'VERSION INFO'  1>&2 ;
            echo \$(java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar --version 1>&2) ;
            echo '/VERSION INFO' 1>&2 ;
            
            java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -nct $GATK_THREADS -R $ESSENTIAL_BWA_REF --dbsnp ${ESSENTIAL_KNOWN_VARIANTS} -I $input -o $output $GATK_FLAGS
        ""","VariantCallHC"
    }
    
}
