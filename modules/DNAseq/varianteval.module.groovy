// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load  PIPELINE_ROOT + "/modules/DNAseq/varianteval.vars.groovy"

VariantEval = {
    doc title: "GATK Base Quality Recalibration",
        desc:  "Recalibrate Base Qualities in BAM files, using GATK",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = QC + '/GATK_varianteval'
    def GATK_FLAGS = " "

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"])

    transform (".vcf.gz") to (".report") {
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php
        exec """
            ${TOOL_ENV} &&

            export TMPDIR=/tmp &&
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi                                       &&

            java -Djava.io.tmpdir=$TMPDIR -jar \${gatk} -T VariantEval -R $GATK_BWA_REF -nt $GATK_THREADS --dbsnp ${GATK_KNOWN_VARIANTS} --eval $input -o $output;
        ""","VariantEval"
    }
    forward input
}

