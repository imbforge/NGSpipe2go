DE_DESeq2_miRNAmature = {
    doc title: "DE_DESeq2_miRNAmature",
        desc:  "Differential expression analysis based only on mature miRNA counts using linear models and DESeq2",
        constraints: "Only simple contrasts in 1-factor design. Include always the intercept. Always gene filtering",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Anke Busch"

    var subdir : ""
    output.dir = DE_DESeq2_miRNAmature_vars.outdir + "/$subdir"

    def DE_DESeq2_miRNAmature_FLAGS =
        (DE_DESeq2_miRNAmature_vars.targets   ? " targets="   + DE_DESeq2_miRNAmature_vars.targets   : "") +
        (DE_DESeq2_miRNAmature_vars.contrasts ? " contrasts=" + DE_DESeq2_miRNAmature_vars.contrasts : "") +
        (DE_DESeq2_miRNAmature_vars.filter    ? " filter="    + DE_DESeq2_miRNAmature_vars.filter    : "") +
        (DE_DESeq2_miRNAmature_vars.prefix    ? " prefix="    + DE_DESeq2_miRNAmature_vars.prefix    : "") +
        (DE_DESeq2_miRNAmature_vars.suffix    ? " suffix="    + DE_DESeq2_miRNAmature_vars.suffix    : "") +
        (DE_DESeq2_miRNAmature_vars.cwd       ? " cwd="       + DE_DESeq2_miRNAmature_vars.cwd    + "/$subdir" : "") +
        (DE_DESeq2_miRNAmature_vars.outdir    ? " out="       + DE_DESeq2_miRNAmature_vars.outdir + "/$subdir" : "") +
        (DE_DESeq2_miRNAmature_vars.genes     ? " gtf="       + DE_DESeq2_miRNAmature_vars.genes     : "") +
        (DE_DESeq2_miRNAmature_vars.pattern   ? " pattern="   + DE_DESeq2_miRNAmature_vars.pattern   : "") +
        (DE_DESeq2_miRNAmature_vars.FDR       ? " FDR="       + DE_DESeq2_miRNAmature_vars.FDR       : "") +
        (DE_DESeq2_miRNAmature_vars.FC        ? " FC="        + DE_DESeq2_miRNAmature_vars.FC        : "") +
        (DE_DESeq2_miRNAmature_vars.extra     ? " "           + DE_DESeq2_miRNAmature_vars.extra     : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The DE_DESeq2 module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if DE_DESeq2 should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("DE_DESeq2.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/DE_DESeq2/DE_DESeq2_miRNAmature.R $DE_DESeq2_miRNAmature_FLAGS
        ""","DE_DESeq2_miRNAmature"
    }
}

