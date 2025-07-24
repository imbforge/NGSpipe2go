readAggrData_bioc = {
    doc title: "readAggrData_bioc",
        desc:  "Read cellranger aggr data from single cell experiment",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = readAggrData_bioc_vars.outdir
    
    // in case of ScaleBio use external result directory
    aggr_data_dir = readAggrData_bioc_vars.seqtype == "ScaleBio" ? EXTERNAL_RESULTDIR + "/" + readAggrData_bioc_vars.aggr_data_dir : readAggrData_bioc_vars.aggr_data_dir

    println "aggr_data_dir: " + aggr_data_dir

    def readAggrData_bioc_FLAGS =
        (readAggrData_bioc_vars.targets            ? " targets="            + readAggrData_bioc_vars.targets            : "") +
        (readAggrData_bioc_vars.outdir             ? " outdir="             + readAggrData_bioc_vars.outdir             : "") +
        (readAggrData_bioc_vars.pipeline_root      ? " pipeline_root="      + readAggrData_bioc_vars.pipeline_root      : "") +
        (readAggrData_bioc_vars.res                ? " res="                + readAggrData_bioc_vars.res                : "") +
        (readAggrData_bioc_vars.seqtype            ? " seqtype="            + readAggrData_bioc_vars.seqtype            : "") +
        (readAggrData_bioc_vars.gtf                ? " gtf="                + readAggrData_bioc_vars.gtf                : "") +
        (readAggrData_bioc_vars.aggr_data_dir      ? " aggr_data_dir="      + aggr_data_dir                                : "") +
        (readAggrData_bioc_vars.run_demux          ? " run_demux="          + readAggrData_bioc_vars.run_demux          : "") +
        (readAggrData_bioc_vars.demux_out          ? " demux_out="          + readAggrData_bioc_vars.demux_out          : "") +
        (readAggrData_bioc_vars.demuxCluster_out   ? " demuxCluster_out="   + readAggrData_bioc_vars.demuxCluster_out   : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The readAggrData_bioc module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if readAggrData_bioc should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    produce("readAggrData_bioc.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null &&

            Rscript ${PIPELINE_ROOT}/tools/sc_readData/readAggrData_bioc.R $readAggrData_bioc_FLAGS
        ""","readAggrData_bioc"
    }
}

