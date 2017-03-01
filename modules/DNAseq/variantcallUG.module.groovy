//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
VariantCallUG = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"

	output.dir = RESULTS
    def GATK_FLAGS = ""
    
    // check if a region limit was provided
    if (GATK_CALL_REGION!=null && GATK_CALL_REGION.length()>0) {
        GATK_FLAGS = " -L " + GATK_CALL_REGION
    } else {
        GATK_FLAGS = ""
    }
    
    transform (".duprm.realigned.recalibrated.bam") to (".UG.vcf.gz") {
    // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php
        exec """
           
            module load jdk/${JAVA_VERSION} &&

            if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                       &&

            java -Djava.io.tmpdir=$TMPDIR -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt $GATK_THREADS -nct $GATK_THREADS -R $GATK_BWA_REF -glm BOTH -I $input -o $output $GATK_FLAGS
        ""","VariantCallUG"
    }
    
}
