load MODULE_FOLDER + "DNAseq/variantcallHC.vars.groovy"

VariantCallHC = {
    doc title: "GATK Variant Calling HC",
        desc:  "Call variants in BAM files using GATK HaplotypeCaller",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Oliver Drechsel"

    output.dir = RESULTS
    def GATK_FLAGS = ""

    // check if a region limit was provided
    if (GATK_CALL_REGION!=null && GATK_CALL_REGION.length()>0) {
        GATK_FLAGS = " -L " + GATK_CALL_REGION
    } else {
        GATK_FLAGS = ""
    }

    transform (".duprm.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        exec """
            module load jdk/${JAVA_VERSION} &&
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi &&

            java -Djava.io.tmpdir=$TMPDIR -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T HaplotypeCaller -nct $GATK_THREADS -R $GATK_BWA_REF --dbsnp $GATK_KNOWN_VARIANTS -I $input -o $output $GATK_FLAGS;
        ""","VariantCallHC"
    }
}

