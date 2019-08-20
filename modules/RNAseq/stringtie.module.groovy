// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/stringtie.vars.groovy"

StringTie = {
    doc title: "STRING_TIE transcript quantification",
        desc:  "Quantifiy transcripts ",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9",
        author: "Nastasja Kreim"

    output.dir = STRINGTIE_OUT

    // create the LOGS/STRING_TIE folder if it doesn't exists
    def F_LOG = new File(LOGS + "/STRING_TIE")
    if(! F_LOG.exists()) {
        F_LOG.mkdirs()
    }

    // star flags
    def STRINGTIE_FLAGS = STRINGTIE_GTF + " " +
                          STRINGTIE_THREADS + " " +
                          STRINGTIE_EXTRA

    // no|yes|reverse
    if (STRINGTIE_STRANDED == "yes") {
        STRINGTIE_FLAGS = "--fr " + STRINGTIE_FLAGS
    }
    else if (STRINGTIE_STRANDED == "reverse") {
        STRINGTIE_FLAGS = "--rf " + STRINGTIE_FLAGS
    }

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

