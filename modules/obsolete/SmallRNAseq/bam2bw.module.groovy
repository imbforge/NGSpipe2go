load MODULE_FOLDER + "SmallRNAseq/bam2bw.vars.groovy"

Bam2bw = {
    doc title: "Bam2bw",
        desc:  "Convert BAM file to bigWig",
        constraints: "none.",
        author: "Sergi Sayols, Antonio Domingues"

    output.dir=TRACKS

    def SAMPLE_NAME = input.split("/")[-1].split("\\.", 2)[0]

    transform(".bam") to (".scaled.bw") {
        exec """
            module load bedtools/${BEDTOOLS_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&
            module load kentUtils/${KENTUTILS_VERSION} &&

            if [ ! -d ${TMP} ]; then
                mkdir -p ${TMP};
            fi &&

            if [ "${NORMALIZATION_TO_NONSTRUCT}" == "yes" ];
            then
                TOTAL_MAPPED=\$(grep $SAMPLE_NAME ${output.dir}/normalization_factors.txt | cut -f 2);
            else
                TOTAL_MAPPED=\$( samtools flagstat $input | head -n1| cut -f1 -d" ");
            fi &&

            CHRSIZES=${TMP}/\$(basename ${input.prefix}).bam2bw.chrsizes &&
            samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&

            SCALE=\$(echo "1000000/\$TOTAL_MAPPED" | bc -l) &&
            genomeCoverageBed -bg -split -scale \${SCALE} -ibam ${input} -g \${CHRSIZES} > ${output.prefix}.bedgraph &&
            bedGraphToBigWig ${output.prefix}.bedgraph \${CHRSIZES} $output &&
            rm \${CHRSIZES} ${output.prefix}.bedgraph
        ""","Bam2bw"
    }
}

