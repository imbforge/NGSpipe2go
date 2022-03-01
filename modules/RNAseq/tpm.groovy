tpm = {
    doc title: "tpm",
        desc:  "Calculation TPMs based on raw counts",
        constraints: "",
        bpipe_version: "",
        author: "Anke Busch"

    output.dir  = tpm_vars.outdir + "/$subdir"
    def TPM_FLAGS =
        (tpm_vars.genesgtf ? " -g " + tpm_vars.genesgtf : "") +
        (tpm_vars.feature  ? " -f " + tpm_vars.feature  : "") +
        (tpm_vars.extra    ? " "    + tpm_vars.extra    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".readcounts.tsv") to (".tpm.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

             Rscript ${PIPELINE_ROOT}/tools/TPMs/TPMs.R -c $input -o $output $TPM_FLAGS
    
        ""","tpm"
    }
}
