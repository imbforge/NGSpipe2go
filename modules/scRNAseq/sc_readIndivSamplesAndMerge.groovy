sc_readIndivSamplesAndMerge = {
    doc title: "sc_readIndivSamplesAndMerge",
        desc:  "Read individual data from single cell experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sivarajan Karunanithi"

    output.dir = sc_readIndivSamplesAndMerge_vars.outdir
    
    def sc_readIndivSamplesAndMerge_FLAGS =
        (sc_readIndivSamplesAndMerge_vars.targets            ? " targets="            + sc_readIndivSamplesAndMerge_vars.targets            : "") +
        (sc_readIndivSamplesAndMerge_vars.outdir             ? " outdir="             + sc_readIndivSamplesAndMerge_vars.outdir             : "") +
        (sc_readIndivSamplesAndMerge_vars.project            ? " project="            + sc_readIndivSamplesAndMerge_vars.project            : "") +
        (sc_readIndivSamplesAndMerge_vars.res                ? " res="                + sc_readIndivSamplesAndMerge_vars.res                : "") +
        (sc_readIndivSamplesAndMerge_vars.org                ? " org="                + sc_readIndivSamplesAndMerge_vars.org                : "") +
        (sc_readIndivSamplesAndMerge_vars.db                 ? " db="                 + sc_readIndivSamplesAndMerge_vars.db                 : "") +
        (sc_readIndivSamplesAndMerge_vars.gtf                ? " gtf="                + sc_readIndivSamplesAndMerge_vars.gtf                : "") +
        (sc_readIndivSamplesAndMerge_vars.mtgenes            ? " mtgenes="            + sc_readIndivSamplesAndMerge_vars.mtgenes            : "") +
        (sc_readIndivSamplesAndMerge_vars.selgenes           ? " selgenes="           + sc_readIndivSamplesAndMerge_vars.selgenes           : "") +
        (sc_readIndivSamplesAndMerge_vars.cellranger_output  ? " cellranger_output="  + sc_readIndivSamplesAndMerge_vars.cellranger_output  : "") +
        (sc_readIndivSamplesAndMerge_vars.run_demux          ? " run_demux="          + sc_readIndivSamplesAndMerge_vars.run_demux          : "") +
        (sc_readIndivSamplesAndMerge_vars.demux_out          ? " demux_out="          + sc_readIndivSamplesAndMerge_vars.demux_out          : "") +
        (sc_readIndivSamplesAndMerge_vars.demuxCluster_out   ? " demuxCluster_out="   + sc_readIndivSamplesAndMerge_vars.demuxCluster_out   : "") +
        (sc_readIndivSamplesAndMerge_vars.colorByFactor      ? " colorByFactor="      + sc_readIndivSamplesAndMerge_vars.colorByFactor      : "") +
        (sc_readIndivSamplesAndMerge_vars.extra              ? " "                    + sc_readIndivSamplesAndMerge_vars.extra              : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_readIndivSamplesAndMerge module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_readIndivSamplesAndMerge should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_readIndivSamplesAndMerge.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_readData/sc_readIndivSamplesAndMerge.R $sc_readIndivSamplesAndMerge_FLAGS
        ""","sc_readIndivSamplesAndMerge"
    }
}

