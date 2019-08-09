// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/bam2bw.vars.groovy"

bam2bw = {
    doc title: "bam2bw",
        desc:  "Convert BAM file to bigWig",
        constraints: "none.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir=TRACKS

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("kentutils", tools["kentutils"]["version"], tools["kentutils"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])

    transform(".bam") to ("_scaled.bw") {
        exec """
            ${TOOL_ENV} &&

            if [ -n "\$SLURM_JOBID" ]; then
                                export TMPDIR=/jobdir/\${SLURM_JOBID};
                        fi &&

            BASEOUTPUT=`basename $output` &&

            CHRSIZES=\${TMPDIR}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
                        samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
                        TOTAL_MAPPED=\$( samtools flagstat $input | head -n5 | tail -n1 | cut -f1 -d" ") &&
                        SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
                        genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} | sortBed -i -  > \${TMPDIR}/\${BASEOUTPUT%.bw}.bedgraph &&
                        bedGraphToBigWig \${TMPDIR}/\${BASEOUTPUT%.bw}.bedgraph \${CHRSIZES} \${TMPDIR}/\${BASEOUTPUT} &&
                        cp \${TMPDIR}/\${BASEOUTPUT} $output
        ""","bam2bw"
    }
}

