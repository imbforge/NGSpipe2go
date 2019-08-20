// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/strandSpecificBW.vars.groovy"

strandBigWig = {
    doc title: "strandBigWig",
        desc:  "strandBigWig wrapper",
        constraints: "strandspecific bigwig for rnaseq",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = STRANDSPECIFICBIGWIG_OUTDIR
    //this might be confusing regarding the reverse and forward
    //setting but according to the deeptools manual it has to be like
    //that. 
     if(ESSENTIAL_STRANDED == "yes") {
        def FORWARD="reverse"
        def REVERSE="forward"
    } else if(ESSENTIAL_STRANDED == "reverse") {
        def FORWARD="forward"
        def REVERSE="reverse"
    }

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
            bamCoverage --numberOfProcessors $STRANDSPECIFICBIGWIG_CORES --filterRNAstrand $FORWARD $STRANDSPECIFICBIGWIG_OTHER -b $input -o \${TMP}/\${base}.fwd.bedgraph &&
            bamCoverage --numberOfProcessors $STRANDSPECIFICBIGWIG_CORES --filterRNAstrand $REVERSE $STRANDSPECIFICBIGWIG_OTHER -b $input -o \${TMP}/\${base}.rev.bedgraph &&

            awk 'BEGIN {OFS="\t"; FS="\t"}{print \$1,\$2, \$3,"-"\$4}' \${TMP}/\${base}.rev.bedgraph > \${TMP}/\${base}.rev.bedgraph_neg &&
            mv \${TMP}/\${base}.rev.bedgraph_neg \${TMP}/\${base}.rev.bedgraph &&
            bedGraphToBigWig \${TMP}/\${base}.fwd.bedgraph $CHRSIZES $output1 &&
            bedGraphToBigWig \${TMP}/\${base}.rev.bedgraph $CHRSIZES $output2
        ""","strandBigWig"
    }
    forward input
}

