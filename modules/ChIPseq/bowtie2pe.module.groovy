// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2pe.vars.groovy"

bowtie2_pe = {
    doc title: "Bowtie PE alignment",
        desc:  "Align paired end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = BOWTIE2_MAPPED

    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /.R1.fastq.gz/).replaceFirst("")


    def BOWTIE2_FLAGS = "-q "  +
                       BOWTIE2_QUALS    + " " + 
                       BOWTIE2_MM_SEED  + " " + 
                       BOWTIE2_INSERT   + " " + 
                       BOWTIE2_THREADS  + " " + 
                       BOWTIE2_EXTRA
    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " + BOWTIE2_SAMTOOLS_THREADS

    def TOOL_ENV = prepare_tool_env("bowtie2", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie2_pe")

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bowtie2 $BOWTIE2_FLAGS $BOWTIE2_REF -1 $input1 -2 $input2 | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output;
        ""","bowtie2_pe"
    }
}



