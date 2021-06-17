STAR_pe_2nd = {
    doc title: "STAR PE alignment",
       desc:  "Align single end and pair end reads",
       constraints: "Only works with compressed input. Set all global vars.",
       author: "Sergi Sayols modified by N.Kreim and Antonio Domingues"

    output.dir = STAR_pe_2nd_vars.outdir

    //before we start we have to define the output filename
    File f = new File(input1)
    def EXP = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    def STAR_FLAGS =
        " --runMode alignReads "        +
        " --genomeLoad NoSharedMemory " +
        " --outStd SAM "                +
        " --outSAMattributes Standard " +
        " --outSJfilterReads Unique "   +
        " --readFilesCommand zcat "     +
        " --outFileNamePrefix "         + STAR_pe_2nd_vars.outdir + "/" + EXP + "." +
        (STAR_pe_2nd_vars.maxram   ? " --limitGenomeGenerateRAM " + STAR_pe_2nd_vars.maxram   : "") +
        (STAR_pe_2nd_vars.bufsize  ? " --limitIObufferSize "      + STAR_pe_2nd_vars.bufsize  : "") +
        (STAR_pe_2nd_vars.ref      ? " --genomeDir "              + STAR_pe_2nd_vars.ref      : "") +
        (STAR_pe_2nd_vars.threads  ? " --runThreadN "             + STAR_pe_2nd_vars.threads  : "") +
        (STAR_pe_2nd_vars.mm       ? " --outFilterMismatchNmax "  + STAR_pe_2nd_vars.mm       : "") +
        (STAR_pe_2nd_vars.multimap ? " --outFilterMultimapNmax "  + STAR_pe_2nd_vars.multimap : "") +
        (STAR_pe_2nd_vars.minintro ? " --alignIntronMin "         + STAR_pe_2nd_vars.minintro : "") +
        (STAR_pe_2nd_vars.overhang ? " --sjdbOverhang "           + STAR_pe_2nd_vars.overhang : "") +
        (STAR_pe_2nd_vars.gtf      ? " --sjdbGTFfile "            + STAR_pe_2nd_vars.gtf      : "") +
        (STAR_pe_2nd_vars.extra    ? " "                          + STAR_pe_2nd_vars.extra    : "")

    // samtools flags
    def SAMTOOLS_VIEW_FLAGS = "-bhSu" +
        (STAR_pe_2nd_vars.filter_sec ? " -F 256" : "")
    def SAMTOOLS_SORT_FLAGS = " -O bam" +
        (STAR_pe_2nd_vars.samtools_threads ? " -@ " + STAR_pe_2nd_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("star", tools["star"]["version"], tools["star"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:"", input:new File(input1.prefix).getName())

    produce(EXP + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            STAR $STAR_FLAGS --outTmpDir \${TMP}/$EXP --readFilesIn $inputs | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/${EXP}_sort - > $output &&
            rm -rf \${TMP}/${EXP}*
        ""","STAR_pe_2nd"
   }
}
