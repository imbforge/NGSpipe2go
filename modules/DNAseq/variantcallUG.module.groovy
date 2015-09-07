//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
VariantCallUG = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	author: "Oliver Drechsel"

	output.dir = RESULTS
    
    transform (".realigned.recalibrated.bam") to (".UG.vcf.gz") {
        
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php
        def GATK_FLAGS = " --output_mode EMIT_ALL_SITES " 
        
        exec """
        
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
        
        java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt $GATK_THREADS -nct $GATK_THREADS -R $ESSENTIAL_BWA_REF -glm BOTH -I $input -o $output
        
        ""","VariantCallUG"
    }
    
}