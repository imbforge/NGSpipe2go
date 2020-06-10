pbc = {
    doc title: "PBC",
        desc:  "PCR Bottleneck Coefficient",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = pbc_vars.outdir

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("pbc")

    transform(".bam") to("_PBC.csv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/ENCODEqc/PBC.R $input && mv ${input.prefix}_PBC.csv $output.dir
        ""","pbc"
    }

    forward input

}
