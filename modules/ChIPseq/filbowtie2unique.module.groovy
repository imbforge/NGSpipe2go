// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/filbowtie2unique.vars.groovy"

filbowtie2unique = {
    doc title: "filter out multimapping reads from bowtie2 out",
        desc:  "filter out multimapping reads from bowtie2 out. output bam file",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = FILBOWTIE2UNIQUE_MAPPED

    def TOOL_ENV = prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("filbowtie2unique")

    transform(".bam") to (".unique.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            samtools view -f 2 $FILBOWTIE2UNIQUE_SAMTOOLS_MAPQ -bhu ${input} | samtools sort $FILBOWTIE2UNIQUE_SAMTOOLS_THREADS -T \${TMP}/\$(basename $output.prefix) -o ${output} -;
        ""","filbowtie2unique"
    }
}

