sc_readAggrData = {
    doc title: "sc_readAggrData",
        desc:  "Read cellranger aggr data from single cell experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_readAggrData_vars.outdir
    
    def sc_readAggrData_FLAGS =
        (sc_readAggrData_vars.targets            ? " targets="            + sc_readAggrData_vars.targets            : "") +
        (sc_readAggrData_vars.outdir             ? " outdir="             + sc_readAggrData_vars.outdir             : "") +
        (sc_readAggrData_vars.project            ? " project="            + sc_readAggrData_vars.project            : "") +
        (sc_readAggrData_vars.res                ? " res="                + sc_readAggrData_vars.res                : "") +
        (sc_readAggrData_vars.org                ? " org="                + sc_readAggrData_vars.org                : "") +
        (sc_readAggrData_vars.db                 ? " db="                 + sc_readAggrData_vars.db                 : "") +
        (sc_readAggrData_vars.gtf                ? " gtf="                + sc_readAggrData_vars.gtf                : "") +
        (sc_readAggrData_vars.mtgenes            ? " mtgenes="            + sc_readAggrData_vars.mtgenes            : "") +
        (sc_readAggrData_vars.selgenes           ? " selgenes="           + sc_readAggrData_vars.selgenes           : "") +
        (sc_readAggrData_vars.cellranger_aggr_id ? " cellranger_aggr_id=" + sc_readAggrData_vars.cellranger_aggr_id : "") +
        (sc_readAggrData_vars.run_demux          ? " run_demux="          + sc_readAggrData_vars.run_demux          : "") +
        (sc_readAggrData_vars.demux_out          ? " demux_out="          + sc_readAggrData_vars.demux_out          : "") +
        (sc_readAggrData_vars.demuxCluster_out   ? " demuxCluster_out="   + sc_readAggrData_vars.demuxCluster_out   : "") +
        (sc_readAggrData_vars.colorByFactor      ? " colorByFactor="      + sc_readAggrData_vars.colorByFactor      : "") +
        (sc_readAggrData_vars.extra              ? " "                    + sc_readAggrData_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_readAggrData module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_readAggrData should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_readAggrData.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_readData/sc_readAggrData.R $sc_readAggrData_FLAGS
        ""","sc_readAggrData"
    }
}

