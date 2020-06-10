trackhub_config = {
    doc title: "Trackhub Configuration",
        desc:  "Generate UCSC track hub configuration file for project tracks",
        constraints: "Currently only processes ChIP-seq & RNA-seq tracks (read from the 'targets.txt' file)",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Martin Oti"

    def TRACKHUB_FLAGS =
        (trackhub_config_vars.project    ? "ESSENTIAL_PROJECT="    + trackhub_config_vars.project    : "") +
        (trackhub_config_vars.db         ? "ESSENTIAL_DB="         + trackhub_config_vars.db         : "") +
        (trackhub_config_vars.chromsizes ? "ESSENTIAL_CHROMSIZES=" + trackhub_config_vars.chromsizes : "") +
        (trackhub_config_vars.stranded   ? "ESSENTIAL_STRANDED="   + trackhub_config_vars.stranded   : "") +
        (trackhub_config_vars.targets    ? "TRACKHUB_TARGETS="     + trackhub_config_vars.targets    : "") +
        (trackhub_config_vars.ftpbase    ? "TRACKHUB_FTPBASE="     + trackhub_config_vars.ftpbase    : "") +
        (trackhub_config_vars.ftpurlbase ? "TRACKHUB_FTPURLBASE="  + trackhub_config_vars.ftpurlbase : "") +
        (trackhub_config_vars.ucsccfg    ? "TRACKHUB_UCSCCFG="     + trackhub_config_vars.ucsccfg    : "") +
        (trackhub_config_vars.peaksdir   ? "TRACKHUB_PEAKSDIR="    + trackhub_config_vars.peaksdir   : "") +
        (trackhub_config_vars.tracksdir  ? "TRACKHUB_TRACKSDIR="   + trackhub_config_vars.tracksdir  : "") +
        (trackhub_config_vars.config     ? "TRACKHUB_CONFIG="      + trackhub_config_vars.config     : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("trackhub_config")

    produce(TRACKHUB_CONFIG) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/trackhub/Configure_Trackhub.R $TRACKHUB_FLAGS
        ""","trackhub_config"
    }
}

