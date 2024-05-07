grn = {
    doc title: "grn",
        desc:  "Identify gene regulatory networks with Pando",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = grn_vars.outdir
    
    def grn_FLAGS =
        (grn_vars.outdir             ? " outdir="             + grn_vars.outdir             : "") +
        (grn_vars.project            ? " project="            + grn_vars.project            : "") +
        (grn_vars.res                ? " res="                + grn_vars.res                : "") +
        (grn_vars.peak_assay         ? " peak_assay="         + grn_vars.peak_assay         : "") +
        (grn_vars.rna_assay          ? " rna_assay="          + grn_vars.rna_assay          : "") +
        (grn_vars.db                 ? " db="                 + grn_vars.db                 : "") +
        (grn_vars.methodModel        ? " methodModel="        + grn_vars.methodModel        : "") +
        (grn_vars.genes2use          ? " genes2use="          + grn_vars.genes2use          : "") +
        (grn_vars.pval_thresh        ? " pval_thresh="        + grn_vars.pval_thresh        : "") +
        (grn_vars.min_genes          ? " min_genes="          + grn_vars.min_genes          : "") +
        (grn_vars.features4graph     ? " features4graph="     + grn_vars.features4graph     : "") +
        (grn_vars.umap_method        ? " umap_method="        + grn_vars.umap_method        : "") +
        (grn_vars.n_neighbors        ? " n_neighbors="        + grn_vars.n_neighbors        : "") +
        (grn_vars.batchCorrection    ? " batchCorrection="    + grn_vars.batchCorrection    : "") +
        (grn_vars.extra              ? " "                    + grn_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The grn module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if grn should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("grn.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_grn/grn.R $grn_FLAGS
        ""","grn"
    }
}

