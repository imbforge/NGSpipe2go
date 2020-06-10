bowtie2 = {
    doc title: "Bowtie2 alignment",
        desc:  "Align single or paired end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = bowtie2_vars.mapped
   
    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    def BOWTIE2_INPUT = (bowtie2_vars.paired   ?  " -1 $input1 -2 $input2.optional" : " -U $input")

    def BOWTIE2_FLAGS = "-q " + 
        (bowtie2_vars.quals      ?  " "     + bowtie2_vars.quals      : "") + 
        (bowtie2_vars.presets    ?  " "     + bowtie2_vars.presets    : "") + 
        (bowtie2_vars.seedMM     ?  " -N "  + bowtie2_vars.seedMM     : "") +
        (bowtie2_vars.seedSubstr ?  " -L "  + bowtie2_vars.seedSubstr : "") +
        (bowtie2_vars.interval   ?  " -i "  + bowtie2_vars.interval   : "") +
        (bowtie2_vars.extend     ?  " -D "  + bowtie2_vars.extend     : "") +
        (bowtie2_vars.sets       ?  " -R "  + bowtie2_vars.sets       : "") +
        (bowtie2_vars.threads    ?  " -p "  + bowtie2_vars.threads    : "") + 
        (bowtie2_vars.extra      ?  " "     + bowtie2_vars.extra      : "") +
        (bowtie2_vars.ref        ?  " -x "  + bowtie2_vars.ref        : "") +
        (bowtie2_vars.paired     ?  " "     + bowtie2_vars.pe_vars    : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie2_vars.samtools_threads ? " -@ " + bowtie2_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie2", tools["bowtie2"]["version"], tools["bowtie2"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie2")

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bowtie2 $BOWTIE2_FLAGS $BOWTIE2_INPUT | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output;
        ""","bowtie2"
    }
}



