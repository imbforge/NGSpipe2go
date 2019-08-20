// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1.vars.groovy"

bowtie_se = {
    doc title: "Bowtie SE alignment",
        desc:  "Align single end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = MAPPED

    def BOWTIE_FLAGS = "-q --sam "  +
                       BOWTIE_QUALS    + " " +
                       BOWTIE_BEST     + " " +
                       BOWTIE_MM_SEED  + " " +
                       BOWTIE_INSERT   + " " +
                       BOWTIE_MAQERR   + " " +
                       BOWTIE_MULTIMAP + " " +
                       BOWTIE_THREADS  + " " +
                       BOWTIE_EXTRA

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " + BOWTIE_SAMTOOLS_THREADS

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie_se")

    transform(".fastq.gz") to (".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix)_bowtie1_sort - > $output
        ""","bowtie_se"
    }
}

