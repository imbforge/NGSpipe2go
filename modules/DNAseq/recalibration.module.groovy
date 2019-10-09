// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load  PIPELINE_ROOT + "/modules/DNAseq/recalibration.vars.groovy"

BaseRecalibration = {
    doc title: "GATK Base Quality Recalibration",
        desc:  "Recalibrate Base Qualities in BAM files, using GATK",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = BaseRecalibration_vars.outdir

    def BaseRecalibrator_FLAGS =
        (BaseRecalibration_vars.known_variants ? " -knownSites " + BaseRecalibration_vars.known_variants : "" ) +
        (BaseRecalibration_vars.threads        ? " -nct "        + BaseRecalibration_vars.threads        : "" ) +
        (BaseRecalibration_vars.bwa_ref        ? " -R "          + BaseRecalibration_vars.bwa_ref        : "" )

    def PrintReads_FLAGS =
        (BaseRecalibration_vars.threads ? " -nct " + BaseRecalibration_vars.threads : "" ) +
        (BaseRecalibration_vars.bwa_ref ? " -R "   + BaseRecalibration_vars.bwa_ref     : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("BaseRecalibration")

    transform (".bam") to (".recalibration.table", ".recalibrated.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${BaseRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T BaseRecalibrator $BaseRecalibrator_FLAGS -I $input -o $output1 &&
            java ${BaseRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T PrintReads $PrintReads_FLAGS -I $input -BQSR $output1 -o $output2
        ""","BaseRecalibration"
    }
    
}
