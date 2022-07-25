// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

pear = {
    doc title: "pear",
        desc:  "Merge overlapping reads from Illumina paired-end sequencing (see Zhang et al (2014) Bioinformatics 30(5):614-620. doi:10.1093/bioinformatics/btt593)",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Frank Rühle"

    output.dir  = pear_vars.outdir

    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")


    // create the log folder if it doesn't exists
    def pear_logdir = new File(pear_vars.logdir)
    if (!pear_logdir.exists()) {
        pear_logdir.mkdirs()
    }
    // create the unassembled folder if it doesn't exists
    def pear_unassembled_dir = new File(pear_vars.outdir + "/unassembled")
    if (!pear_unassembled_dir.exists()) {
        pear_unassembled_dir.mkdirs()
    }


    def PEAR_FLAGS =
        (pear_vars.pvalue ? " -p " + pear_vars.pvalue : "") +
        (pear_vars.minoverlap  ? " -v " + pear_vars.minoverlap : "") +
        (pear_vars.quality_threshold  ? " -q " + pear_vars.quality_threshold : "") +
        (pear_vars.min_trim_length  ? " -t " + pear_vars.min_trim_length : "") +
        (pear_vars.max_uncalled_base  ? " -u " + pear_vars.max_uncalled_base : "") +
        (pear_vars.test_method  ? " -g " + pear_vars.test_method : "") +
        (pear_vars.score_method  ? " -s " + pear_vars.score_method : "") +
        (pear_vars.phred_base  ? " -b " + pear_vars.phred_base : "") +
        (pear_vars.cap  ? " -c " + pear_vars.cap : "") +
        (pear_vars.threads  ? " -j " + pear_vars.threads  : "") +
        (pear_vars.extra  ? " "             + pear_vars.extra       : "")

    def TOOL_ENV = prepare_tool_env("pear", tools["pear"]["version"], tools["pear"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // run the chunk
    produce(OUTPUTFILE + ".assembled.fastq.gz"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            pear $PEAR_FLAGS -f $input1 -r $input2 -o \${TMP}/\$(basename ${OUTPUTFILE}) 1> $pear_logdir/\$(basename ${OUTPUTFILE}).log &&
            
            gzip \${TMP}/\$(basename ${OUTPUTFILE}).assembled.fastq && 
            mv \${TMP}/\$(basename ${OUTPUTFILE}).assembled.fastq.gz $output.dir/ &&

            gzip \${TMP}/\$(basename ${OUTPUTFILE}).[ud]*.fastq && 
            mv \${TMP}/\$(basename ${OUTPUTFILE}).[ud]*.fastq.gz $pear_unassembled_dir 

        ""","pear"
    }
}


