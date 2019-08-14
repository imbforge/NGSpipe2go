// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load  PIPELINE_ROOT + "/modules/DNAseq/variantfuseHC.vars.groovy"

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

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])

    transform (".duprm.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        exec """
            ${TOOL_ENV} &&

            export TMPDIR=/tmp &&
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi                                       &&

            java -Djava.io.tmpdir=$TMPDIR -jar \${gatk} -T HaplotypeCaller -nct $GATK_THREADS -R $GATK_BWA_REF --dbsnp $GATK_KNOWN_VARIANTS -I $input -o $output;
        ""","VariantFuseHC"
    }
}

