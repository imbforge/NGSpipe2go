StringTie = {
    doc title: "STRING_TIE transcript quantification",
        desc:  "Quantifiy transcripts ",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9",
        author: "Nastasja Kreim"

    output.dir = StringTie_vars.outdir

    // create the LOGS/STRING_TIE folder if it doesn't exists
    def F_LOG = new File(LOGS + "/STRING_TIE")
    if(! F_LOG.exists()) {
        F_LOG.mkdirs()
    }

    // star flags
    def STRINGTIE_FLAGS =
        (StringTie_gtf      ? " -G " + stringtie_gtf     : "" ) +
        (StringTie_threads  ? " -p " + stringtie_threads : "" ) +
        (StringTie_extra    ? " "    + stringtie_extra   : "" ) 

    // no|yes|reverse
    if (StringTie_stranded == "yes") {
        STRINGTIE_FLAGS = "--fr " + STRINGTIE_FLAGS
    }
    else if (StringTie_stranded == "reverse") {
        STRINGTIE_FLAGS = "--rf " + STRINGTIE_FLAGS
    }
    // NOTE: else? what if stranded == "no"?

    def TOOL_ENV = prepare_tool_env("stringtie", tools["stringtie"]["version"], tools["stringtie"]["runenv"])
    def PREAMBLE = get_preamble("StringTie")

    // code chunk
    transform(".bam") to("_stringtie.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            base=\$(basename $input) &&
            base=\${base%.bam} &&
            echo \$base &&
            echo $STRINGTIE_FLAGS &&
            mkdir $output.dir/\$base &&
            stringtie $input $STRINGTIE_FLAGS -A $output.dir/\${base}/\${base}_gene_abund.tab -o $output.dir/\$base/\${base}_stringtie.gtf &&
            touch $output
        ""","StringTie"
      }
}

