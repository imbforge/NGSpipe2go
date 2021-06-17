geneBodyCov = {
    doc title: "geneBodyCoverage",
        desc: "Calculate the RNA-seq coverage over gene body. Useful to check the 5' or 3' coverage bias",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = geneBodyCov_vars.outdir
    def GENEBODYCOV_FLAGS =
        (geneBodyCov_vars.format ? " -f " + geneBodyCov_vars.format : "" ) +
        (geneBodyCov_vars.bed    ? " -r " + geneBodyCov_vars.bed    : "" ) +
        (geneBodyCov_vars.extra  ? " "    + geneBodyCov_vars.extra  : "" ) 

    def TOOL_ENV = prepare_tool_env("rseqc", tools["rseqc"]["version"], tools["rseqc"]["runenv"])
    def PREAMBLE = get_preamble(module:"geneBodyCoverage", branch:branch, branch_outdir:"")

    // run the chunk
    transform(".bam") to (".geneBodyCoverage.curves.png", ".geneBodyCoverage.r", ".geneBodyCoverage.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            geneBody_coverage.py -i $input -o ${output3.prefix.prefix} $GENEBODYCOV_FLAGS
        ""","geneBodyCov"
    }
    forward input
}
