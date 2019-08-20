// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/scRNAseq/subread.vars.groovy"

subread_count = {
    doc title: "subread_count",
        desc:  "Counting reads in features with feature-count out of the subread package, modified to allow umitools specific counting",
        constraints: """Default: strand specific counting.""",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel, Nastasja Kreim"

    output.dir  = SUBREAD_OUTDIR
    def SUBREAD_FLAGS = "--donotsort" +  " " + 
                        SUBREAD_CORES    + " " +
                        SUBREAD_GENESGTF + " " +
                        SUBREAD_EXTRA    + " " +
                        "-R BAM" + " " + 
                        "--Rpath " + output.dir 

    if(SUBREAD_PAIRED == "yes") {
        SUBREAD_FLAGS = "-p " + SUBREAD_FLAGS
    }

    // no|yes|reverse
    if(SUBREAD_STRANDED == "no") {
        SUBREAD_FLAGS = "-s 0 " + SUBREAD_FLAGS
    }
    else if (SUBREAD_STRANDED == "yes") {
        SUBREAD_FLAGS = "-s 1 " + SUBREAD_FLAGS
    }
    else {
        SUBREAD_FLAGS = "-s 2 " + SUBREAD_FLAGS
    }

    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("subread_count")

    // run the chunk
    transform(".bam") to (".featureCounts.bam", ".raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            featureCounts $SUBREAD_FLAGS -o $output2 $input 2> ${output.prefix}_subreadlog.stderr &&
            base=`basename $input` &&
            samtools sort ${output.dir}/\${base}.featureCounts.bam > $output1 &&
            rm ${output.dir}/\${base}.featureCounts.bam;
        ""","subread_count"
    }
}
