trackhub = {
    doc title: "Trackhub",
        desc:  "Generate UCSC track hub to display project tracks",
        constraints: "Uses configuration file, which should have been generated earlier",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Martin Oti"

    def TRACKHUB_FLAGS =
        (trackhub_vars.config ? "TRACKHUB_CONFIG=" + trackhub_vars.config : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"]) + " && " +
                   prepare_tool_env("kentutils", tools["kentutils"]["version"], tools["kentutils"]["runenv"])
    def PREAMBLE = get_preamble("trackhub")

    transform(".yaml") to (".done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/trackhub/Make_Trackhub.R $TRACKHUB_FLAGS
        ""","trackhub"
    }
}

