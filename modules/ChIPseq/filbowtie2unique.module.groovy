// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.vars.groovy"

filbowtie2unique = {
    doc title: "filter out multimapping reads from bowtie2 out",
        desc:  "filter out multimapping reads from bowtie2 out. output bam file",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = filbowtie2unique_vars.mapped

    def FILBOWTIE2_FLAGS = (filbowtie2unique_vars.paired ? " -f 2 -q $filbowtie2unique_vars.samtools_mapq_pe" : " -F 4 -q $filbowtie2unique_vars.samtools_mapq_se")

    def TOOL_ENV = prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("filbowtie2unique")

    transform(".bam") to (".unique.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            samtools view $FILBOWTIE2_FLAGS -bhu ${input} | samtools sort -@ $filbowtie2unique_vars.samtools_threads -T \${TMP}/\$(basename $output.prefix) -o ${output} -
        ""","filbowtie2unique"
    }
}

