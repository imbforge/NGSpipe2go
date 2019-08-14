// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/subread.vars.groovy"

subread_count = {
    doc title: "subread_count_se",
        desc:  "Counting reads in features with feature-count out of the subread package",
        constraints: """Default: strand specific counting.""",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir  = SUBREAD_OUTDIR
    def SUBREAD_FLAGS = "--donotsort" +  " " + 
                        SUBREAD_CORES    + " " +
                        SUBREAD_GENESGTF + " " +
                        SUBREAD_EXTRA    + " "

    if(SUBREAD_PAIRED == "yes") {
        SUBREAD_FLAGS = "-p " + SUBREAD_FLAGS
    }

    // no|yes|reverse
    if(SUBREAD_STRANDED == "no") {
        SUBREAD_FLAGS = "-s 0 " + SUBREAD_FLAGS
    }
    else if (SUBREAD_STRANDED == "yes") {
        SUBREAD_FLAGS = "-s 1 " + SUBREAD_FLAGS
    }
    else {
        SUBREAD_FLAGS = "-s 2 " + SUBREAD_FLAGS
    }

    def TOOL_ENV = prepare_tool_env("subread", tools["subread"]["version"], tools["subread"]["runenv"])

    // run the chunk
    transform(".bam") to (".raw_readcounts.tsv") {
        exec """
            ${TOOL_ENV} &&
    
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi &&
            base=`basename $input` &&
            if [[ "$SUBREAD_PAIRED" == "yes" ]]; 
            then
                echo "We are resorting and doing the repair\n" &&
                repair -i $input $SUBREAD_CORES -o \${TMPDIR}/\${base} &&
                featureCounts $SUBREAD_FLAGS -o $output \${TMPDIR}/\${base} 2> ${output.prefix}_subreadlog.stderr;
            else
                featureCounts $SUBREAD_FLAGS -o $output $input 2> ${output.prefix}_subreadlog.stderr;
            fi
        ""","subread_count"
    }
}
