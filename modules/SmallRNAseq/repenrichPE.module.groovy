RepEnrichPE = {
    doc title: "Quantification of Transposons",
        desc:  "Quantifies transposon expression/targeting using RepEnrich. First does unique mapping to the genome keeping the multimapped reads in a fastq file (--max). Both the uniquely mapped and multimapping reads are then used for repEnrich. Intermediate results are removed.",
        constraints: "BedTools <= 2.23.0",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Antonio Domingues"

    def EXP = input1.split("/")[-1].replaceAll("_R1.fastq.gz", "")
    output.dir = REPENR_OUT_DIR + "/" + EXP
    def INPUT_TMP1 = output.dir + "/" + EXP + "_R1.fastq"
    def INPUT_TMP2 = output.dir + "/" + EXP + "_R2.fastq"

    produce (
            EXP + ".bt.log",  
            EXP + "_fraction_counts.txt") {

    exec """
        module load RepEnrich/0.1 &&
        module load bowtie/${BOWTIE_VERSION} &&
        module load bedtools/2.23.0 &&
        module load samtools/${SAMTOOLS_VERSION} &&

        MULTI=$output.dir/${EXP}".multimap.fastq" 
        UNIQ=$output.dir/${EXP}".bam"

        echo "Uncompress files:" &&

        zcat $input1 > $INPUT_TMP1 &&
        zcat $input2 > $INPUT_TMP2 &&

        echo "Mapping:" &&


        bowtie -p $REPENRICH_CORES $BOWTIE_RE_OTHER --max $output.dir/${EXP}".multimap.fastq"  $ESSENTIAL_BOWTIE_REF -1 $INPUT_TMP1 -2 $INPUT_TMP2 2> $output1 | samtools view -bhSu - | samtools sort -@ $REPENRICH_CORES -o $output.dir/${EXP}".bam" -T ${EXP} - &&
        samtools index $output.dir/${EXP}".bam" &&

        echo "Repeat counting:" &&
        
        python RepEnrich.py ${REPEAT_MASKER} $output.dir ${EXP} ${REPEAT_REF} $output.dir/${EXP}.multimap_1.fastq --fastqfile2 $output.dir/${EXP}.multimap_2.fastq $output.dir/${EXP}".bam" --cpus ${REPENRICH_CORES} --is_bed ${REPENRICH_BED} --pairedend TRUE &&

        rm $INPUT_TMP1 $INPUT_TMP2

        ""","RepEnrichPE"
    }
}
