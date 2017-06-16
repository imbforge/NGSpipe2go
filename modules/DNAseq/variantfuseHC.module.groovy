//rule for task GATK recalibration from DNAseq pipeline
//desc: Base quality recalibration provided by GATK
VariantFuseHC = {
    
    doc title: "GATK Base Quality Recalibration",
	desc:  "Recalibrate Base Qualities in BAM files, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"

	output.dir = RESULTS
    def GATK_FLAGS = GATK_REFCONF   + " " +
                     GATK_INDEXTYPE + " " +
                     GATK_INDEXPARM + " " +
                     GATK_EXTRA
    
    transform (".duprm.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        
        exec """
            module load jdk/${JAVA_VERSION} &&

            if [ -n "\$SLURM_JOBID" ]; then
				export TMPDIR=/jobdir/\${SLURM_JOBID} &&
                       		mkdir ${TMPDIR};
			fi                                       &&

            java -Djava.io.tmpdir=$TMPDIR -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T HaplotypeCaller -nct $GATK_THREADS -R $GATK_BWA_REF --dbsnp $GATK_KNOWN_VARIANTS -I $input -o $output &&
            rm -rf ${TMPDIR};
        ""","VariantFuseHC"
    }
}
