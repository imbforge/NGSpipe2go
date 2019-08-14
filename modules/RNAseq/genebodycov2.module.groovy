// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/genebodycov2.vars.groovy"

geneBodyCov2 = {
    doc title: "geneBodyCoverage2",
        desc:  """Calculate the RNA-seq coverage over gene body. 
            Useful to check the 5' or 3' coverage bias""",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.9",
        author: "Sergi Sayols"

    output.dir = GENEBODYCOV2_OUTDIR.replaceFirst("outdir=", "")
    def GENEBODYCOV2_FLAGS = GENEBODYCOV2_GTF      + " " +
                             GENEBODYCOV2_PAIRED   + " " +
                             GENEBODYCOV2_STRANDED + " " +
                             GENEBODYCOV2_OUTDIR   + " " +
                             GENEBODYCOV2_THREADS

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    // run the chunk
    transform(".bam") to ("_geneBodyCov.png") {
        exec """
            ${TOOL_ENV} &&

            if [[ ! -e "$output.dir" ]]; then
                mkdir -p "$output.dir";
            fi &&

            Rscript ${PIPELINE_ROOT}/tools/geneBodyCov/geneBodyCov.R bam=$input $GENEBODYCOV2_FLAGS
        ""","geneBodyCov2"
    }
    forward input
}
