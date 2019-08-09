// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

pbc = {
    doc title: "PBC",
        desc:  "PCR Bottleneck Coefficient",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = QC + "/pbc"

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    transform(".bam") to("_PBC.csv") {
        exec """
            ${TOOL_ENV} &&

            Rscript ${PIPELINE_ROOT}/tools/ENCODEqc/PBC.R $input && mv ${input.prefix}_PBC.csv $output.dir
        ""","pbc"
    }

    forward input

}
