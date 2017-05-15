RepEnrich = {
    doc title: "Quantification of Transposons",
        desc:  "Quantifies transposon expression/targeting using RepEnrich. First does unique mapping to the genome keeping the multimapped reads in a fastq file (--max). Both the uniquely mapped and multimapping reads are then used for repEnrich. Intermediate results are removed.",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Antonio Domingues"

    def SAMPLE_NAME = input.split("/")[-1].replaceAll(".deduped_barcoded.trimmed.fastq.gz", "")
    output.dir = REPENR_OUT_DIR + "/" + SAMPLE_NAME

    def BOWTIE_FLAGS = " -q --sam -t"  +
                    " "   + BOWTIE_QUALS    +
                    " "   + BOWTIE_BEST     +
                    " -p " + Integer.toString(BOWTIE_THREADS) +
                    " -m " + Integer.toString(BOWTIE_MULTIMAP) 

    produce (
            SAMPLE_NAME + ".bt.log",  
            SAMPLE_NAME + "_fraction_counts.txt") {

    exec """
        echo 'VERSION INFO'  1>&2 &&
        ${TOOL_BOWTIE}/bowtie --version 1>&2 &&
        echo '/VERSION INFO' 1>&2 &&

        MULTI=${SAMPLE_NAME}".multimap.fastq" 
        UNIQ=${SAMPLE_NAME}".bam"

        zcat $input | ${TOOL_BOWTIE}/bowtie $BOWTIE_FLAGS --max ${MULTI} $BOWTIE_REF - 2> $output1 | ${TOOL_SAMTOOLS}/samtools view -bhSu - | ${TOOL_SAMTOOLS}/samtools sort -@ $BOWTIE_THREADS - -o ${UNIQ} &&
        ${TOOL_SAMTOOLS}/samtools index ${UNIQ} &&

        python ${TOOL_REPENRICH}/RepEnrich.py ${REPEAT_MASKER} $output.dir ${SAMPLE_NAME} ${REPEAT_REF} ${MULTI} ${UNIQ} --cpus ${BOWTIE_THREADS} --is_bed ${REPENRICH_BED} &&

        rm ${MULTI} ${UNIQ}

        ""","RepEnrich"
    }
}
