load MODULE_FOLDER + "SmallRNAseq/repenrich.vars.groovy"

RepEnrich = {
    doc title: "Quantification of Transposons",
        desc:  "Quantifies transposon expression/targeting using RepEnrich. First does unique mapping to the genome keeping the multimapped reads in a fastq file (--max). Both the uniquely mapped and multimapping reads are then used for repEnrich. Intermediate results are removed.",
        constraints: "BedTools <= 2.23.0",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Antonio Domingues"

    def SAMPLE_NAME = input.split("/")[-1].replaceAll(".fastq.gz", "")
    output.dir = REPENR_OUT_DIR + "/" + SAMPLE_NAME

    def REPENRICH_FLAGS = " -q --sam -t"  +
                    " "   + REPENRICH_QUALS    +
                    " "   + REPENRICH_BEST     +
                    " -p " + Integer.toString(REPENRICH_ESSENTIAL_THREADS) +
                    " -m " + Integer.toString(REPENRICH_MULTIMAP) 

    produce (
            SAMPLE_NAME + ".bt.log",  
            SAMPLE_NAME + "_fraction_counts.txt") {

    exec """
        module load RepEnrich/${REPENRICH_VERSION} &&
        module load bowtie/${BOWTIE_VERSION} &&
        module load bedtools/2.23.0 &&
        module load samtools/${SAMTOOLS_VERSION} &&

        MULTI=${SAMPLE_NAME}".multimap.fastq" 
        UNIQ=${SAMPLE_NAME}".bam"

        zcat $input | bowtie $REPENRICH_FLAGS --max $output.dir/${MULTI} $REPENRICH_REF - 2> $output1 | samtools view -bhSu - | samtools sort -@ $REPENRICH_ESSENTIAL_THREADS - -o $output.dir/${UNIQ} &&
        samtools index $output.dir/${UNIQ} &&

        RepEnrich.py ${REPEAT_MASKER} $output.dir ${SAMPLE_NAME} ${REPEAT_REF} $output.dir/${MULTI} $output.dir/${UNIQ} --cpus ${REPENRICH_THREADS} --is_bed ${REPENRICH_BED} &&

        rm $output.dir/${MULTI} $output.dir/${UNIQ}

        ""","RepEnrich"
    }
}
