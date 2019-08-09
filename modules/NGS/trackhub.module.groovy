// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/trackhub.vars.groovy"

trackhub = {
    doc title: "Trackhub",
        desc:  "Generate UCSC track hub to display project tracks",
        constraints: "Uses configuration file, which should have been generated earlier",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Martin Oti"

    def trackhub_FLAGS = "TRACKHUB_CONFIG=" + TRACKHUB_CONFIG

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"]) + " && " +
                   prepare_tool_env("kentutils", tools["kentutils"]["version"], tools["kentutils"]["runenv"])

    transform(".yaml") to (".done") {
        exec """
            ${TOOL_ENV} &&

            Rscript ${PIPELINE_ROOT}/tools/trackhub/Make_Trackhub.R $trackhub_FLAGS
        ""","trackhub"
    }
}

