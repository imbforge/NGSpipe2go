// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/trackhub_config.vars.groovy"

trackhub_config = {
    doc title: "Trackhub Configuration",
        desc:  "Generate UCSC track hub configuration file for project tracks",
        constraints: "Currently only processes ChIP-seq & RNA-seq tracks (read from the 'targets.txt' file)",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Martin Oti"

    def trackhub_FLAGS = "ESSENTIAL_PROJECT="    + ESSENTIAL_PROJECT    + " " +
                         "ESSENTIAL_DB="         + ESSENTIAL_DB         + " " +
                         "ESSENTIAL_CHROMSIZES=" + ESSENTIAL_CHROMSIZES + " " +
                         "ESSENTIAL_STRANDED="   + ESSENTIAL_STRANDED   + " " +
                         "TRACKHUB_TARGETS="     + TRACKHUB_TARGETS     + " " +
                         "TRACKHUB_FTPBASE="     + TRACKHUB_FTPBASE     + " " +
                         "TRACKHUB_FTPURLBASE="  + TRACKHUB_FTPURLBASE  + " " +
                         "TRACKHUB_UCSCCFG="     + TRACKHUB_UCSCCFG     + " " +
                         "TRACKHUB_PEAKSDIR="    + TRACKHUB_PEAKSDIR    + " " +
                         "TRACKHUB_TRACKSDIR="   + TRACKHUB_TRACKSDIR   + " " +
                         "TRACKHUB_CONFIG="      + TRACKHUB_CONFIG

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("trackhub_config")

    produce(TRACKHUB_CONFIG) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/trackhub/Configure_Trackhub.R $trackhub_FLAGS
        ""","trackhub_config"
    }
}

