//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
bowtie_se = {
    doc title: "Bowtie SE alignment",
        desc:  "Align single end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = MAPPED

    def BOWTIE_FLAGS = "-q --sam "  +
                       BOWTIE_QUALS    + " " + 
                       BOWTIE_BEST     + " " + 
                       BOWTIE_MM_SEED  + " " + 
                       BOWTIE_INSERT   + " " + 
                       BOWTIE_MAQERR   + " " + 
                       BOWTIE_MULTIMAP + " " + 
                       BOWTIE_THREADS  + " " + 
                       BOWTIE_EXTRA

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " + BOWTIE_SAMTOOLS_THREADS

    transform(".fastq.gz") to (".bam") {
        exec """
            module load bowtie/${BOWTIE_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&

            zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T $TMP/\$(basename $output.prefix)_bowtie1_sort - > $output
        ""","bowtie_se"
    }
}

