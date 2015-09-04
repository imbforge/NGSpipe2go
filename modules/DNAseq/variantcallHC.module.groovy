//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
VariantCallHC = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	author: "Oliver Drechsel"

	output.dir = RESULTS
    
    transform (".realigned.recalibrated.bam") to (".HC.vcf.gz") {
        
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        def GATK_FLAGS = "--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000"
        
        
        exec """
        
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
        java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -nct $GATK_THREADS -R $ESSENTIAL_BWA_REF --dbsnp ${ESSENTIAL_KNOWN_VARIANTS} -I $input -o $output
        
        ""","VariantCallHC"
    }
    
}