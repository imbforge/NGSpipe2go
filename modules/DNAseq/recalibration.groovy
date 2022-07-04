BaseRecalibration = {
    doc title: "GATK Base Quality Recalibration",
        desc:  "Recalibrate Base Qualities in BAM files, using GATK",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel, modified by Frank Rühle"

    output.dir = BaseRecalibration_vars.outdir

    // create folder for recalibration table if it doesn't exists
    def BQSR_STATSDIR = new File( BaseRecalibration_vars.statsdir)
    if (!BQSR_STATSDIR.exists()) {
        BQSR_STATSDIR.mkdirs()
    }
    def SAMPLENAME = new File(input)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_TABLE = SAMPLENAME_BASE.replace(".bam", ".recalibration.table") 
    def SAMPLENAME_BASE_PLOT  = SAMPLENAME_BASE.replace(".bam", ".analyzeCovariates.pdf") 

    def BaseRecalibrator_FLAGS =
        (BaseRecalibration_vars.known_variants ? " --known-sites " + BaseRecalibration_vars.known_variants : "" ) +
        (BaseRecalibration_vars.bwa_ref        ? " -R "            + BaseRecalibration_vars.bwa_ref        : "" ) + 
        (BaseRecalibration_vars.extra          ? " "               + BaseRecalibration_vars.extra          : "" )

    def ApplyBQSR_FLAGS =
        (BaseRecalibration_vars.bwa_ref ? " -R "   + BaseRecalibration_vars.bwa_ref     : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"]) + " && " +
                   prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform (".bam") to (".recalibrated.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk BaseRecalibrator $BaseRecalibrator_FLAGS -I $input -O ${BQSR_STATSDIR}/${SAMPLENAME_BASE_TABLE} &&
            gatk AnalyzeCovariates -bqsr ${BQSR_STATSDIR}/${SAMPLENAME_BASE_TABLE} -plots ${BQSR_STATSDIR}/${SAMPLENAME_BASE_PLOT} &&
            gatk ApplyBQSR $ApplyBQSR_FLAGS -I $input --bqsr-recal-file ${BQSR_STATSDIR}/${SAMPLENAME_BASE_TABLE} -O $output

        ""","BaseRecalibration"
    }

}
