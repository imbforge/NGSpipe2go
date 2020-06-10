bowtie2_pe = {
    doc title: "Bowtie2 PE alignment",
        desc:  "Align paired end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = bowtie2_pe_vars.mapped

    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /.R1.fastq.gz/).replaceFirst("")

    def BOWTIE2_FLAGS = "-q " +
        (bowtie2_pe_vars.quals   ?  " "     + bowtie2_pe_vars.quals   : "") + 
        (bowtie2_pe_vars.mm_seed ?  " "     + bowtie2_pe_vars.mm_seed : "") + 
        (bowtie2_pe_vars.insert  ?  " -L "  + bowtie2_pe_vars.insert  : "") +
        (bowtie2_pe_vars.threads ?  " -p "  + bowtie2_pe_vars.threads : "") + 
        (bowtie2_pe_vars.extra   ?  " "     + bowtie2_pe_vars.extra   : "") +
        (bowtie2_pe_vars.ref     ?  " -x "  + bowtie2_pe_vars.ref     : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie2_pe_vars.samtools_threads ? " -@ " + bowtie2_pe_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie2", tools["bowtie2"]["version"], tools["bowtie2"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie2_pe")

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bowtie2 $BOWTIE2_FLAGS -1 $input1 -2 $input2 | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output;
        ""","bowtie2_pe"
    }
}



