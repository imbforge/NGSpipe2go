// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/base_recalibration.vars.groovy"

BaseRecalibration = {
    doc title: "GATK BaseRecalibrator",
        desc: "Recalibrate Base Qualities in BAM files, using GATK.",
        constraints: "GATK version >= 3.5",
        author: "Antonio Domingues"

    output.dir = OUTDIR_STAR2ND

    def JAVA_FLAGS = "-Xmx" + RECAL_MAXMEM
    def GATK_FLAGS  = " -R " + GATK_REF +
                      " -nct " + GATK_THREADS

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])

    transform (".bam") to (".recalibration.table", ".recalibrated.bam"){

        exec """
            ${TOOL_ENV} &&

            echo 'VERSION INFO'  1>&2 &&
            echo \$(java -jar \${gatk} --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            java -jar \${gatk} -T BaseRecalibrator -I $input -o $output1 -knownSites $VCF_REF  $GATK_FLAGS &&

            java $JAVA_FLAGS -jar \${gatk} -T PrintReads -I $input -BQSR $output1 -o $output2 $GATK_FLAGS
        ""","BaseRecalibration"
    }
}
