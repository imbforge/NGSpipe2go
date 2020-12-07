bam2bw = {
    doc title: "bam2bw",
        desc:  "Convert BAM file to bigWig",
        constraints: "none.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir=bam2bw_vars.outdir

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("kentutils", tools["kentutils"]["version"], tools["kentutils"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bam2bw")

    transform(".bam") to ("_scaled.bw") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            BASEOUTPUT=`basename $output` &&
            CHRSIZES=\${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
            samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
            TOTAL_MAPPED=\$( samtools flagstat $input | head -n5 | tail -n1 | cut -f1 -d" ") &&
            SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
            genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} | sortBed -i -  > \${TMP}/\${BASEOUTPUT%.bw}.bedgraph &&
            bedGraphToBigWig \${TMP}/\${BASEOUTPUT%.bw}.bedgraph \${CHRSIZES} \${TMP}/\${BASEOUTPUT} &&
            cp -f \${TMP}/\${BASEOUTPUT} $output
        ""","bam2bw"
    }
}

