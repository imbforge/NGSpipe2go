sc_bioc_readAggrData = {
    doc title: "sc_bioc_readAggrData",
        desc:  "Read cellranger aggr data from single cell experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = sc_bioc_readAggrData_vars.outdir
    
    // in case of ScaleBio use external result directory
    aggr_data_dir = sc_bioc_readAggrData_vars.seqtype == "ScaleBio" ? EXTERNAL_RESULTDIR + "/" + sc_bioc_readAggrData_vars.aggr_data_dir : sc_bioc_readAggrData_vars.aggr_data_dir

    println "aggr_data_dir: " + aggr_data_dir

    def SC_BIOC_READAGGRDATA_FLAGS =
        (sc_bioc_readAggrData_vars.targets            ? " targets="            + sc_bioc_readAggrData_vars.targets            : "") +
        (sc_bioc_readAggrData_vars.outdir             ? " outdir="             + sc_bioc_readAggrData_vars.outdir             : "") +
        (sc_bioc_readAggrData_vars.pipeline_root      ? " pipeline_root="      + sc_bioc_readAggrData_vars.pipeline_root      : "") +
        (sc_bioc_readAggrData_vars.res                ? " res="                + sc_bioc_readAggrData_vars.res                : "") +
        (sc_bioc_readAggrData_vars.seqtype            ? " seqtype="            + sc_bioc_readAggrData_vars.seqtype            : "") +
        (sc_bioc_readAggrData_vars.gtf                ? " gtf="                + sc_bioc_readAggrData_vars.gtf                : "") +
        (sc_bioc_readAggrData_vars.aggr_data_dir      ? " aggr_data_dir="      + aggr_data_dir                                : "") +
        (sc_bioc_readAggrData_vars.run_demux          ? " run_demux="          + sc_bioc_readAggrData_vars.run_demux          : "") +
        (sc_bioc_readAggrData_vars.demux_out          ? " demux_out="          + sc_bioc_readAggrData_vars.demux_out          : "") +
        (sc_bioc_readAggrData_vars.demuxCluster_out   ? " demuxCluster_out="   + sc_bioc_readAggrData_vars.demuxCluster_out   : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The sc_bioc_readAggrData module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if sc_bioc_readAggrData should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("sc_bioc_readAggrData.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_readData/sc_bioc_readAggrData.R $SC_BIOC_READAGGRDATA_FLAGS
        ""","sc_bioc_readAggrData"
    }
}

