DNAaccess = {
    doc title: "DNAaccess",
        desc:  "Process the DNA accessibility assay by performing latent semantic indexing",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = DNAaccess_vars.outdir
    
    def DNAaccess_FLAGS =
        (DNAaccess_vars.outdir             ? " outdir="             + DNAaccess_vars.outdir             : "") +
        (DNAaccess_vars.project            ? " project="            + DNAaccess_vars.project            : "") +
        (DNAaccess_vars.res                ? " res="                + DNAaccess_vars.res                : "") +
        (DNAaccess_vars.featureCutoff      ? " featureCutoff="      + DNAaccess_vars.featureCutoff      : "") +
        (DNAaccess_vars.skipFirstLSIcomp   ? " skipFirstLSIcomp="   + DNAaccess_vars.skipFirstLSIcomp   : "") +
        (DNAaccess_vars.extra              ? " "                    + DNAaccess_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The DNAaccess module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if DNAaccess should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("DNAaccess.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_DNAaccess/DNAaccess.R $DNAaccess_FLAGS
        ""","DNAaccess"
    }
}

