CTannoMarker = {
    doc title: "CTannoMarker",
        desc:  "Annotate cell types by marker genes with scType",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sivarajan Karunanithi, Frank Rühle"

    output.dir = CTannoMarker_vars.outdir
    
    def CTannoMarker_FLAGS =
        (CTannoMarker_vars.outdir              ? " outdir="              + CTannoMarker_vars.outdir              : "") +
        (CTannoMarker_vars.project             ? " project="             + CTannoMarker_vars.project             : "") +
        (CTannoMarker_vars.res                 ? " res="                 + CTannoMarker_vars.res                 : "") +
        (CTannoMarker_vars.dbfile              ? " dbfile="              + CTannoMarker_vars.dbfile              : "") +
        (CTannoMarker_vars.tissue              ? " tissue="              + CTannoMarker_vars.tissue              : "") +        
        (CTannoMarker_vars.extra               ? " "                     + CTannoMarker_vars.extra               : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The CTannoMarker module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if CTannoMarker should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("CTannoMarker.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/CTanno/CTannoMarker.R $CTannoMarker_FLAGS
        ""","CTannoMarker"
    }
}

