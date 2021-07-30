geneBodyCov2 = {
    doc title: "geneBodyCoverage2",
        desc:  """Calculate the RNA-seq coverage over gene body. 
            Useful to check the 5' or 3' coverage bias""",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.9",
        author: "Sergi Sayols"

    output.dir = geneBodyCov2_vars.outdir
    def GENEBODYCOV2_FLAGS =
        (geneBodyCov2_vars.gtf      ? " gtf="      + geneBodyCov2_vars.gtf      : "" ) +
        (geneBodyCov2_vars.paired   ? " paired="   + geneBodyCov2_vars.paired   : "" ) +
        (geneBodyCov2_vars.stranded ? " stranded=" + geneBodyCov2_vars.stranded : "" ) +
        (geneBodyCov2_vars.outdir   ? " outdir="   + geneBodyCov2_vars.outdir   : "" ) +
        (geneBodyCov2_vars.threads  ? " threads="  + geneBodyCov2_vars.threads  : "" ) 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    transform(".bam") to ("_geneBodyCov.png") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [[ ! -e "$output.dir" ]]; then
                mkdir -p "$output.dir";
            fi &&

            Rscript ${PIPELINE_ROOT}/tools/geneBodyCov/geneBodyCov.R bam=$input $GENEBODYCOV2_FLAGS
        ""","geneBodyCov2"
    }
    forward input
}
