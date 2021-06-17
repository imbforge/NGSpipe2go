BaseRecalibration = {
    doc title: "GATK BaseRecalibrator",
        desc: "Recalibrate Base Qualities in BAM files, using GATK.",
        constraints: "GATK version >= 3.5",
        author: "Antonio Domingues"

    output.dir = BaseRecalibration_vars.outdir

    def BaseRecalibrator_FLAGS =
        (BaseRecalibration_vars.vcf_ref    ? " -knownSites " + BaseRecalibration_vars.vcf_ref    : "" ) +
        (BaseRecalibration_vars.threads    ? " -nct "        + BaseRecalibration_vars.threads    : "" ) +
        (BaseRecalibration_vars.genome_ref ? " -R "          + BaseRecalibration_vars.genome_ref : "" )

    def PrintReads_FLAGS =
        (BaseRecalibration_vars.threads    ? " -nct " + BaseRecalibration_vars.threads    : "" ) +
        (BaseRecalibration_vars.genome_ref ? " -R "   + BaseRecalibration_vars.genome_ref : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(module:"BaseRecalibration", branch:branch, branch_outdir:"")

    transform (".bam") to (".recalibration.table", ".recalibrated.bam"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${BaseRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T BaseRecalibrator $BaseRecalibrator_FLAGS -I $input -o $output1 &&
            java ${BaseRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T PrintReads $PrintReads_FLAGS -I $input -BQSR $output1 -o $output2
        ""","BaseRecalibration"
    }
}
