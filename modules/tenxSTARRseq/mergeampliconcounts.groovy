MergeAmpliconCounts = {
    doc title: "MergeAmpliconCounts",
        desc:  "Merge 10X library PCR-amplified STARR mRNA counts with endogenous mRNA counts from original 10X sc/snRNA-seq library.",
        constraints: "If Bioconductor SingleCellExperiment object already exists with 10X sc/snRNA-seq counts, the new counts will be added to it. Otherwise they will be added to the CellRanger counts and put in a new SCE object.",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Martin Oti"

    def SCEFILE = new File(MergeAmpliconCounts_vars.sce)
    output.dir = SCEFILE.getParentFile()

    def MERGEAMPLICONCOUNTS_FLAGS =
        (MergeAmpliconCounts_vars.starrgtf  ? " starrgtf="  + MergeAmpliconCounts_vars.starrgtf  : "" ) +
        (MergeAmpliconCounts_vars.umicounts ? " umicounts=" + MergeAmpliconCounts_vars.umicounts : "" ) +
        (MergeAmpliconCounts_vars.aggr      ? " aggr="      + MergeAmpliconCounts_vars.aggr      : "" ) +
        (MergeAmpliconCounts_vars.prefix    ? " prefix="    + MergeAmpliconCounts_vars.prefix    : "" ) +
        (MergeAmpliconCounts_vars.sce       ? " sce="       + MergeAmpliconCounts_vars.sce       : "" )

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("MergeAmpliconCounts")

    // run the chunk
    produce("merge_amplicon_counts.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [[ ! -e "$output.dir" ]]; then
                mkdir -p "$output.dir";
            fi &&

            Rscript ${PIPELINE_ROOT}/tools/MergeSTARRCounts/MergeAmpliconCounts.R $MERGEAMPLICONCOUNTS_FLAGS
            
            touch $output
            
        ""","MergeAmpliconCounts"
    }
    forward inputs
}
