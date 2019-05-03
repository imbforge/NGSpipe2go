//rule for task STAR_se from catalog RNAseq, version 1
//desc: Align single end reads
StringTie = {
    doc title: "STRING_TIE transcript quantification",
        desc:  "Quantifiy transcripts ",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9",
        author: "Nastasja Kreim"

    output.dir = STRINGTIE_OUT

    // create the LOGS/STRING_TIE folder if it doesn't exists
    F_LOG = new File(LOGS + "/STRING_TIE")
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

    // code chunk
    transform(".bam") to("_stringtie.done") {
      exec """
        module load stringtie/${STRINGTIE_VERSION} && 
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
