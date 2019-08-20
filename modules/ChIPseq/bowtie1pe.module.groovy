// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1pe.vars.groovy"

bowtie_pe = {
    doc title: "Bowtie PE alignment",
        desc:  "Align paired end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols modified for paired end by Nastasja Kreim"

    output.dir = BOWTIE_MAPPED
    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /.R1.fastq.gz/).replaceFirst("")

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
    def PREAMBLE = get_preamble("bowtie_pe")

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

      bowtie $BOWTIE_FLAGS $BOWTIE_REF -1 <(zcat $input1) -2 <(zcat $input2) | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output;
    ""","bowtie_pe"
  }
}

