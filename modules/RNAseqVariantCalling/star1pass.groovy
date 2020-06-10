STAR_pe = {
    doc title: "STAR PE alignment",
       desc:  "Align single end and pair end reads",
       constraints: "Only works with compressed input. Set all global vars.",
       author: "Sergi Sayols modified by N.Kreim and Antonio Domingues"

    output.dir = STAR_pe_vars.outdir

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
        " --outFileNamePrefix "         + STAR_pe_vars.outdir + "/" + EXP + "." +
        (STAR_pe_vars.maxram   ? " --limitGenomeGenerateRAM " + STAR_pe_vars.maxram   : "") +
        (STAR_pe_vars.bufsize  ? " --limitIObufferSize "      + STAR_pe_vars.bufsize  : "") +
        (STAR_pe_vars.ref      ? " --genomeDir "              + STAR_pe_vars.ref      : "") +
        (STAR_pe_vars.threads  ? " --runThreadN "             + STAR_pe_vars.threads  : "") +
        (STAR_pe_vars.mm       ? " --outFilterMismatchNmax "  + STAR_pe_vars.mm       : "") +
        (STAR_pe_vars.multimap ? " --outFilterMultimapNmax "  + STAR_pe_vars.multimap : "") +
        (STAR_pe_vars.minintro ? " --alignIntronMin "         + STAR_pe_vars.minintro : "") +
        (STAR_pe_vars.overhang ? " --sjdbOverhang "           + STAR_pe_vars.overhang : "") +
        (STAR_pe_vars.gtf      ? " --sjdbGTFfile "            + STAR_pe_vars.gtf      : "") +
        (STAR_pe_vars.extra    ? " "                          + STAR_pe_vars.extra    : "")

    def TOOL_ENV = prepare_tool_env("star", tools["star"]["version"], tools["star"]["runenv"])
    def PREAMBLE = get_preamble("STAR_pe")

    produce(EXP + ".SJ.out.tab") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            STAR $STAR_FLAGS --outTmpDir \${TMP}/$EXP --readFilesIn $inputs > /dev/null &&
            rm -rf \${TMP}/${EXP}*
        ""","STAR_pe"
    }
}
