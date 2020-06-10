bowtie_pe = {
    doc title: "Bowtie PE alignment",
        desc:  "Align paired end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols modified for paired end by Nastasja Kreim"

    output.dir = bowtie_pe_vars.mapped
    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /.R1.fastq.gz/).replaceFirst("")

    def BOWTIE_FLAGS = " -q --sam"  +
        (bowtie_pe_vars.best    ? " --best --strata --tryhard --chunkmbs 256 " : "") +
        (bowtie_pe_vars.quals   ? " "    + bowtie_pe_vars.quals    : "") +
        (bowtie_pe_vars.mm_seed ? " -n " + bowtie_pe_vars.mm_seed  : "") +
        (bowtie_pe_vars.insert  ? " -l " + bowtie_pe_vars.insert   : "") +
        (bowtie_pe_vars.maqerr  ? " -e " + bowtie_pe_vars.maqerr   : "") +
        (bowtie_pe_vars.multimap_mode && bowtie_pe_vars.multimap ?
          (bowtie_pe_vars.multimap_mode == "discard" ? " -m " : " -M ") + bowtie_pe_vars.multimap : "") +
        (bowtie_pe_vars.threads ? " -p " + bowtie_pe_vars.threads  : "") +
        (bowtie_pe_vars.extra   ? " "    + bowtie_pe_vars.extra    : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie_pe_vars.samtools_threads ? " -@ " + bowtie_pe_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie_pe")

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bowtie $BOWTIE_FLAGS $bowtie_pe_vars.ref -1 <(zcat $input1) -2 <(zcat $input2) | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output;
    ""","bowtie_pe"
  }
}

