motifActivity = {
    doc title: "motifActivity",
        desc:  "Determine differential activity scores between for each cluster or cell type, respectively, compared to all other cluster and celltypes",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = motifActivity_vars.outdir
    
    def motifActivity_FLAGS =
        (motifActivity_vars.outdir             ? " outdir="             + motifActivity_vars.outdir             : "") +
        (motifActivity_vars.project            ? " project="            + motifActivity_vars.project            : "") +
        (motifActivity_vars.res                ? " res="                + motifActivity_vars.res                : "") +
        (motifActivity_vars.db                 ? " db="                 + motifActivity_vars.db                 : "") +
        (motifActivity_vars.clusterVar         ? " clusterVar="         + motifActivity_vars.clusterVar         : "") +
        (motifActivity_vars.CTannoSelected     ? " CTannoSelected="     + motifActivity_vars.CTannoSelected     : "") +
        (motifActivity_vars.motif2plot         ? " motif2plot="         + motifActivity_vars.motif2plot         : "") +
        (motifActivity_vars.extra              ? " "                    + motifActivity_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The motifActivity module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if motifActivity should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("motifActivity.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_motifs/motifActivity.R $motifActivity_FLAGS
        ""","motifActivity"
    }
}

