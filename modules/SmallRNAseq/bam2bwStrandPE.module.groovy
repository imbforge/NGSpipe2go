load MODULE_FOLDER + "SmallRNAseq/bam2bwStrandPE.vars.groovy"
Bam2bwStrandPE = {
    doc title: "Bam2bw",
        desc:  "Convert BAM file to bigWig for each strand. Normalization factor is the total number of reads. First separates reads. Source and inspiration: https://www.biostars.org/p/92935/",
        constraints: "none.",
        author: "Antonio Domingues"

    def EXP = input.split("/")[-1].replaceAll(".bam", "")
    output.dir=BAMCOVSTRANDSPE_OUTDIR
    //this might be confusing regarding the reverse and forward
    //setting but according to the deeptools manual it has to be like
    //that. 
        if(ESSENTIAL_STRANDED == "yes") {
        def FORWARD="antisense"
        def REVERSE="sense"
    } else if(ESSENTIAL_STRANDED == "reverse") {
        def FORWARD="sense"
        def REVERSE="antisense"
    }

    transform(".bam") to (".scaled.sense.bw", ".scaled.antisense.bw")  {
        exec """
            module load bedtools/${BEDTOOLS_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&
            module load kentUtils/${KENTUTILS_VERSION} &&

            if [ ! -d ${TMP} ]; then
                mkdir -p ${TMP};
            fi &&

            CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
            samtools idxstats ${input} | cut -f1-2 > ${CHRSIZES} &&

            TOTAL_MAPPED=\$( samtools flagstat $input | head -n1| cut -f1 -d" ") &&
            SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&


            samtools view -b -f 128 -F 16 $input > ${TMP}/${EXP}.sense1.bam &&
            samtools index ${TMP}/${EXP}.sense1.bam &&

            samtools view -b -f 80 $input > ${TMP}/${EXP}.sense2.bam &&
            samtools index ${TMP}/${EXP}.sense2.bam &&

            samtools merge -f ${TMP}/${EXP}.sense.bam ${TMP}/${EXP}.sense1.bam ${TMP}/${EXP}.sense2.bam &&
            samtools index ${TMP}/${EXP}.sense.bam &&

            samtools view -b -f 144 $input > ${TMP}/${EXP}.antisense1.bam &&
            samtools index ${TMP}/${EXP}.antisense1.bam &&

            samtools view -b -f 64 -F 16 $input > ${TMP}/${EXP}.antisense2.bam &&
            samtools index ${TMP}/${EXP}.antisense2.bam &&

            samtools merge -f ${TMP}/${EXP}.antisense.bam ${TMP}/${EXP}.antisense1.bam ${TMP}/${EXP}.antisense2.bam &&
            samtools index ${TMP}/${EXP}.antisense.bam &&

            genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${TMP}/${EXP}.${FORWARD}.bam -g \${CHRSIZES} > ${output1.prefix}.bedgraph &&
            bedGraphToBigWig ${output1.prefix}.bedgraph \${CHRSIZES} ${output1} &&

            genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${TMP}/${EXP}.${REVERSE}.bam -g \${CHRSIZES} > ${output2.prefix}.bedgraph &&
            bedGraphToBigWig ${output2.prefix}.bedgraph \${CHRSIZES} ${output2} &&

            rm ${CHRSIZES} ${output1.prefix}.bedgraph ${output2.prefix}.bedgraph ${TMP}/${EXP}.*.bam
        ""","Bam2bwStrandPE"
    }
}
