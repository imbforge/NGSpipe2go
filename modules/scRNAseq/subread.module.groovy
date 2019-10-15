// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/scRNAseq/subread.vars.groovy"

subread_count = {
    doc title: "subread_count",
        desc:  "Counting reads in features with feature-count out of the subread package, modified to allow umitools specific counting",
        constraints: """Default: strand specific counting.""",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel, Nastasja Kreim"

    output.dir  = subread_count_vars.outdir
    def SUBREAD_FLAGS =
        "--donotsort " +
        "-R BAM "      + 
        "--Rpath "     + subread_count_vars.outdir +
        (subread_count_vars.threads  ? " -T " + subread_count_vars.threads  : "") +
        (subread_count_vars.genesgtf ? " -a " + subread_count_vars.genesgtf : "") +
        (subread_count_vars.paired   ? " -p "                               : "") +
        (subread_count_vars.extra    ? " "    + subread_count_vars.extra    : "") +
        (subread_count_vars.stranded == "no" ? " -s0 " : (subread2rnatypes_vars.stranded == "yes" ? " -s1 " : " -s2 "))

    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("subread_count")

    // run the chunk
    transform(".bam") to (".featureCounts.bam", ".raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            featureCounts $SUBREAD_FLAGS --tmpDir \${TMP}/featureCounts_${input} -o $output2 $input 2> ${output.prefix}_subreadlog.stderr &&
            base=`basename $input` &&
            samtools sort ${output.dir}/\${base}.featureCounts.bam > $output1 &&
            rm ${output.dir}/\${base}.featureCounts.bam
        ""","subread_count"
    }
}
