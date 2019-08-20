// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastq_quality_filter.vars.groovy"

FastQQualityFilter = {
    doc title: "Remove sequences",
        desc:  "filter reads containing low-quality (Phred score below 20) bases in order to facilitate the PCR duplicates removal.",
        constraints: "Only supports compressed FASTQ files",
        author: "Antonio Domingues, Anke Busch"

    output.dir = FASTQ_QUALITY_FILTER_OUTDIR

    // create the log folder if it doesn't exists
    def FASTQ_QUALITY_FILTER_LOGDIR = new File( FASTQ_QUALITY_FILTER_OUTDIR + "/logs")
    if (!FASTQ_QUALITY_FILTER_LOGDIR.exists()) {
        FASTQ_QUALITY_FILTER_LOGDIR.mkdirs()
    }

    def FASTQ_QUALITY_FILTER_FLAGS= " -q " + MIN_QUAL +
                                    " -p " + MIN_PERCENT +
                                    " -Q " + QUAL_FORMAT +
                                    FASTQ_QUALITY_FILTER_OTHER

    def TOOL_ENV = prepare_tool_env("fastx", tools["fastx"]["version"], tools["fastx"]["runenv"])
    def PREAMBLE = get_preamble("FastQQualityFilter")

    transform(".fastq.gz") to (".highQ.fastq.gz") {
        def SAMPLENAME = input.prefix.prefix    
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            zcat $input | fastq_quality_filter $FASTQ_QUALITY_FILTER_FLAGS -o $output 2>&1 >> ${FASTQ_QUALITY_FILTER_LOGDIR}/\${SAMPLENAME_BASE}.fastq_quality_filter.log
        ""","FastQQualityFilter"
    }
}
