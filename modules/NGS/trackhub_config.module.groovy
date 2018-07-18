trackhub_config = {
    doc title: "Trackhub Configuration",
        desc:  "Generate UCSC track hub configuration file for project tracks",
        constraints: "Currently only processes ChIP-seq & RNA-seq tracks (read from the 'targets.txt' file)",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Martin Oti"

    def trackhub_FLAGS = "ESSENTIAL_PROJECT=" + ESSENTIAL_PROJECT + " " +
        "ESSENTIAL_DB=" + ESSENTIAL_DB + " " +
        "ESSENTIAL_CHROMSIZES=" + ESSENTIAL_CHROMSIZES + " " +
        "ESSENTIAL_STRANDED=" + ESSENTIAL_STRANDED + " " +
        "TRACKHUB_TARGETS=" + TRACKHUB_TARGETS + " " +
        "TRACKHUB_FTPBASE=" + TRACKHUB_FTPBASE + " " +
        "TRACKHUB_FTPURLBASE=" + TRACKHUB_FTPURLBASE + " " +
        "TRACKHUB_UCSCCFG=" + TRACKHUB_UCSCCFG + " " +
        "TRACKHUB_PEAKSDIR=" + TRACKHUB_PEAKSDIR + " " +
        "TRACKHUB_TRACKSDIR=" + TRACKHUB_TRACKSDIR + " " +
        "TRACKHUB_CONFIG=" + TRACKHUB_CONFIG

    produce(ESSENTIAL_PROJECT + "/trackhub.yaml") {
        exec """
           module load R/${R_VERSION} &&
           Rscript ${TOOL_TRACKHUB}/Configure_Trackhub.R $trackhub_FLAGS;
        ""","trackhub_config"
    }
}

