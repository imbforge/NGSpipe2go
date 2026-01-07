kmeans_bioc = {
    doc title: "kmeans_bioc",
        desc:  "Apply k-means clustering on singlecellexperiment object",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = kmeans_bioc_vars.outdir
    
    def kmeans_bioc_FLAGS =
        (kmeans_bioc_vars.outdir             ? " outdir="             + kmeans_bioc_vars.outdir             : "") +
        (kmeans_bioc_vars.seqtype            ? " seqtype="            + kmeans_bioc_vars.seqtype            : "") +
        (kmeans_bioc_vars.pipeline_root      ? " pipeline_root="      + kmeans_bioc_vars.pipeline_root      : "") +
        (kmeans_bioc_vars.res                ? " res="                + kmeans_bioc_vars.res                : "") +
        (kmeans_bioc_vars.data2clust         ? " data2clust="         + kmeans_bioc_vars.data2clust         : "") +
        (kmeans_bioc_vars.kmeans             ? " kmeans="             + kmeans_bioc_vars.kmeans             : "") +
        (kmeans_bioc_vars.annocat_plot       ? " annocat_plot="       + kmeans_bioc_vars.annocat_plot       : "") +
        (kmeans_bioc_vars.annocat_plot2      ? " annocat_plot2="      + kmeans_bioc_vars.annocat_plot2      : "") +
        (kmeans_bioc_vars.plot_pointsize     ? " plot_pointsize="     + kmeans_bioc_vars.plot_pointsize     : "") +
        (kmeans_bioc_vars.plot_pointalpha    ? " plot_pointalpha="    + kmeans_bioc_vars.plot_pointalpha    : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The kmeans_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if kmeans_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("kmeans_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_clust/kmeans_bioc.R $kmeans_bioc_FLAGS
        ""","kmeans_bioc"
    }
}

