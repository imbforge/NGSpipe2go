load MODULE_FOLDER + "SmallRNAseq/bamCoverageStrandsRPKM.vars.groovy"

BamCoverageStrandsRPKM = {
    doc title: "BamCoverage",
        desc:  "Convert BAM file to bigWig for each strand. Normalization factor is the total number of reads.",
        constraints: "Requires deepTools >= 2.2",
        author: "Antonio Domingues"

    output.dir=BAMCOVSTRANDS_OUTDIR

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

    transform(".bam") to (".RPKM.fwd.bw", ".RPKM.rev.bw") {
         exec """
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
            module load deepTools/$DEEPTOOLS_VERSION &&
            module load samtools/$SAMTOOLS_VERSION &&

            TOTAL_MAPPED=\$( samtools flagstat $input | head -n1 | cut -f1 -d" ") &&
            SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&

            bamCoverage --bam $input -o $output1 \
                $BAMCOVERAGE_OTHER \
                --numberOfProcessors $ESSENTIAL_CORES \
                --filterRNAstrand forward &&

            bamCoverage --bam $input -o $output2 \
                $BAMCOVERAGE_OTHER \
                --numberOfProcessors $ESSENTIAL_CORES \
                --filterRNAstrand reverse 
                ""","BamCoverageStrandsRPKM"
    }
}

