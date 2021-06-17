bowtie1 = {
    doc title: "Bowtie1 alignment",
        desc:  "Align single end or paired end reads",
        constraints: "Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Sergi Sayols, Nastasja Kreim, Frank Rühle"

    output.dir = bowtie1_vars.mapped

    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    //def BOWTIE_INPUT = (bowtie1_vars.paired   ?  " -1 <(zcat $input1) -2 <(zcat $input2.optional)" : " <(zcat $input)")
    def BOWTIE_INPUT = (bowtie1_vars.paired   ?  " -1 $input1 -2 $input2.optional" : " $input") // bowtie1 does allow compressed fastq files

    def BOWTIE_FLAGS = "-q --sam "  +
        (bowtie1_vars.quals   ? " "    + bowtie1_vars.quals    : "") +
        (bowtie1_vars.best    ? " --best --strata --tryhard --chunkmbs 256" : "") +
        (bowtie1_vars.mm_seed ? " -n " + bowtie1_vars.mm_seed  : "") +
        (bowtie1_vars.insert  ? " -l " + bowtie1_vars.insert   : "") +
        (bowtie1_vars.maqerr  ? " -e " + bowtie1_vars.maqerr   : "") +
        (bowtie1_vars.multimap_mode && bowtie1_vars.multimap ?
        (bowtie1_vars.multimap_mode == "discard" ? " -m " : " -M ") + bowtie1_vars.multimap : "") +
        (bowtie1_vars.threads ? " -p " + bowtie1_vars.threads  : "") +
        (bowtie1_vars.extra   ? " "    + bowtie1_vars.extra    : "") +
        (bowtie1_vars.paired  ? " "    + bowtie1_vars.pe_vars    : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie1_vars.samtools_threads ? " -@ " + bowtie1_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:"", input:new File(input1.prefix).getName())

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bowtie $BOWTIE_FLAGS ${bowtie1_vars.ref} $BOWTIE_INPUT | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix)_bowtie1_sort - > $output
        ""","bowtie1"
    }
}

