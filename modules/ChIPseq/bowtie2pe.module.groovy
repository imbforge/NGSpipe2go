//desc: Align paired end reads
bowtie2_pe = {
    doc title: "Bowtie PE alignment",
        desc:  "Align paired end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = BOWTIE2_MAPPED
    
    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")


    def BOWTIE2_FLAGS = "-q "  +
                       BOWTIE2_QUALS    + " " + 
                       BOWTIE2_MM_SEED  + " " + 
                       BOWTIE2_INSERT   + " " + 
                       BOWTIE2_THREADS  + " " + 
                       BOWTIE2_EXTRA
    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " + BOWTIE2_SAMTOOLS_THREADS

    produce(OUTPUTFILE + ".bam") {
        exec """
            module load bowtie2/${BOWTIE2_VERSION}   &&
            module load samtools/${SAMTOOLS_VERSION} &&
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi                                       &&
            
            bowtie2 $BOWTIE2_FLAGS $BOWTIE2_REF -1 $input1 -2 $input2 | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T $TMPDIR/\$(basename $output.prefix) - > $output;
        ""","bowtie2_pe"
    }
}



