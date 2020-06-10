strandBigWig = {
    doc title: "strandBigWig",
        desc:  "strandBigWig wrapper",
        constraints: "strandspecific bigwig for rnaseq",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = strandBigWig_vars.outdir

    //this might be confusing regarding the reverse and forward
    //setting but according to the deeptools manual it has to be like this
    if(strandBigWig_vars.stranded == "yes") {
        def FORWARD="reverse"
        def REVERSE="forward"
    } else if(strandBigWig_vars.stranded == "reverse") {
        def FORWARD="forward"
        def REVERSE="reverse"
    }

    // if you use > v3 modify the vars to --normalizeUsing RPKM since the API of deeptools changed
    def STRANDBIGWIG_FLAGS =
        (strandBigWig_vars.threads               ? " --numberOfProcessors " + strandBigWig_vars.threads       : "") +
        (strandBigWig_vars.smoothLength          ? " --smoothLength "       + strandBigWig_vars.smoothLength  : "") +
        (strandBigWig_vars.binSize               ? " --binSize "            + strandBigWig_vars.binSize       : "") +
        (strandBigWig_vars.normalizeUsingRPKM    ? " --normalizeUsingRPKM"                                    : "") +
        (strandBigWig_vars.skipNonCoveredRegions ? " --skipNonCoveredRegions"                                 : "") +
        (strandBigWig_vars.extra                 ? " "                      + strandBigWig_vars.extra         : "")

    def TOOL_ENV = prepare_tool_env("deeptools", tools["deeptools"]["version"], tools["deeptools"]["runenv"]) + " && " +
                   prepare_tool_env("kentutils", tools["kentutils"]["version"], tools["kentutils"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("strandBigWig")

    transform(".bam") to(".fwd.bw", ".rev.bw") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            base=\$(basename $input.prefix) &&
            echo \$base &&
            CHRSIZES=\${TMP}/\${base}.bam2bw.chrsizes &&
            samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
            bamCoverage --filterRNAstrand $FORWARD $STRANDBIGWIG_FLAGS -b $input -o \${TMP}/\${base}.fwd.bedgraph &&
            bamCoverage --filterRNAstrand $REVERSE $STRANDBIGWIG_FLAGS -b $input -o \${TMP}/\${base}.rev.bedgraph &&

            awk 'BEGIN {OFS="\t"; FS="\t"}{print \$1,\$2, \$3,"-"\$4}' \${TMP}/\${base}.rev.bedgraph > \${TMP}/\${base}.rev.bedgraph_neg &&
            mv \${TMP}/\${base}.rev.bedgraph_neg \${TMP}/\${base}.rev.bedgraph &&
            bedGraphToBigWig \${TMP}/\${base}.fwd.bedgraph $CHRSIZES $output1 &&
            bedGraphToBigWig \${TMP}/\${base}.rev.bedgraph $CHRSIZES $output2
        ""","strandBigWig"
    }
    forward input
}

