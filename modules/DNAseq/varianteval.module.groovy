//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
VariantEval = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	author: "Oliver Drechsel"

	output.dir = RESULTS
    
    transform (".vcf.gz") to (".report") {
        
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php
        def GATK_FLAGS = " "
        
        exec """
        
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
        
        java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T VariantEval -R $ESSENTIAL_BWA_REF -nt $GATK_THREADS --eval $input -o $output

        ""","VariantEval"
    }
    forward $input
    
}