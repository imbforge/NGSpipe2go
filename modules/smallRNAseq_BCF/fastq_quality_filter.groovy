FastQQualityFilter = {
    doc title: "Remove sequences",
        desc:  "filter reads containing low-quality (Phred score below 20) bases in order to facilitate the PCR duplicates removal.",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch"

    output.dir = FastQQualityFilter_vars.outdir

    // create the log folder if it doesn't exists
    def FASTQ_QUALITY_FILTER_LOGDIR = new File(FastQQualityFilter_vars.logdir)
    if (!FASTQ_QUALITY_FILTER_LOGDIR.exists()) {
        FASTQ_QUALITY_FILTER_LOGDIR.mkdirs()
    }

    def FASTQ_QUALITY_FILTER_FLAGS=
        (FastQQualityFilter_vars.min_qual    ? " -q " + FastQQualityFilter_vars.min_qual    : "") +
        (FastQQualityFilter_vars.min_percent ? " -p " + FastQQualityFilter_vars.min_percent : "") +
        (FastQQualityFilter_vars.qual_format ? " -Q " + FastQQualityFilter_vars.qual_format : "") +
        (FastQQualityFilter_vars.extra       ? " "    + FastQQualityFilter_vars.extra       : "")

    def TOOL_ENV = prepare_tool_env("fastx", tools["fastx"]["version"], tools["fastx"]["runenv"])
    def PREAMBLE = get_preamble(module:"FastQQualityFilter", branch:branch, branch_outdir:"")

    transform(".fastq.gz") to (".highQ.fastq.gz") {
        def SAMPLENAME = input.prefix.prefix    
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            zcat $input | fastq_quality_filter $FASTQ_QUALITY_FILTER_FLAGS -o $output 2>&1 >> ${FastQQualityFilter_vars.logdir}/\${SAMPLENAME_BASE}.fastq_quality_filter.log
        ""","FastQQualityFilter"
    }
}
