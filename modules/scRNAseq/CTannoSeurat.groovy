CTannoSeurat = {
    doc title: "CTannoSeurat",
        desc:  "Annotate cell types by transferring cell labels from an existing reference dataset with Seurat",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = CTannoSeurat_vars.outdir
    
    def CTannoSeurat_FLAGS =
        (CTannoSeurat_vars.outdir              ? " outdir="              + CTannoSeurat_vars.outdir              : "") +
        (CTannoSeurat_vars.project             ? " project="             + CTannoSeurat_vars.project             : "") +
        (CTannoSeurat_vars.res                 ? " res="                 + CTannoSeurat_vars.res                 : "") +
        (CTannoSeurat_vars.pathRefDataset      ? " pathRefDataset="      + CTannoSeurat_vars.pathRefDataset      : "") +
        (CTannoSeurat_vars.columnNameCelltypes ? " columnNameCelltypes=" + CTannoSeurat_vars.columnNameCelltypes : "") +        
        (CTannoSeurat_vars.extra               ? " "                     + CTannoSeurat_vars.extra               : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The CTannoSeurat module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if CTannoSeurat should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("CTannoSeurat.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/CTanno/CTannoSeurat.R $CTannoSeurat_FLAGS
        ""","CTannoSeurat"
    }
}

