// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/bamcoverage.vars.groovy"

bamCoverage = {
    doc title: "bamCoverage",
        desc:  "bamCoverage wrapper",
        constraints: "normalised bigwig track for RNA/ChipSeq PE data",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

        output.dir = BAMCOVERAGE_OUTDIR
        def BAMCOVERAGE_FLAGS = BAMCOVERAGE_CORES + " " + BAMCOVERAGE_OTHER
        if( BAMCOVERAGE_FRAGMENTS == "yes") {
            BAMCOVERAGE_FLAGS = BAMCOVERAGE_FLAGS + " --extendReads"
        }

    def TOOL_ENV = prepare_tool_env("deeptools", tools["deeptools"]["version"], tools["deeptools"]["runenv"])

    transform(".bam") to(".bw") {
        exec """
            ${TOOL_ENV} &&
    
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi;

            bamCoverage $BAMCOVERAGE_FLAGS --bam $input -o ${output};
        ""","bamCoverage"
    }
}

