load MODULE_FOLDER + "DNAseq/bwa.vars.groovy"

BWA_pe = {
    doc title: "BWA PE alignment",
    desc:  "Align paired end reads",
    constraints: "Only works with compressed input. Set all global vars. Samtools version >= 1.2.",
    bpipe_version: "tested with bpipe 0.9.8.7",
    author: "Oliver Drechsel, Anke Busch"

    output.dir = MAPPED

    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")

    def BWA_FLAGS = "-M " +
                    BWA_THREADS + " " +
                    BWA_EXTRA

    def SAMTOOLS_VIEW_FLAGS = "-bhSu"
    def SAMTOOLS_SORT_FLAGS = SAMTOOLS_THREADS

    produce(OUTPUTFILE + ".bam") {
        exec """
            module load bwa/${BWA_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} && 

            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi;

            SAMPLE_NAME=\$(basename $output.prefix) &&
            PLATFORM="genomics" &&

            bwa mem $BWA_FLAGS -R \"@RG\\tID:\${SAMPLE_NAME}\\tSM:\${SAMPLE_NAME}\\tPL:illumina\\tLB:\${SAMPLE_NAME}\\tPU:\${PLATFORM}\" $BWA_REF $input1 $input2 | samtools view ${SAMTOOLS_VIEW_FLAGS} - | samtools sort ${SAMTOOLS_SORT_FLAGS} -T \${TMPDIR}/\${SAMPLE_NAME} -  > ${output} &&

            samtools flagstat ${output} 1>&2
            ""","BWA_pe"
    }
}


BWA_se = {
    doc title: "BWA SE alignment",
    desc:  "Align paired end reads",
    constraints: "Only works with compressed input. Set all global vars. Samtools version >= 1.2",
    bpipe_version: "tested with bpipe 0.9.8.7",
    author: "Oliver Drechsel"

    output.dir = MAPPED
    def BWA_FLAGS = "-M " +
                    BWA_THREADS + " " +
                    BWA_EXTRA

    def SAMTOOLS_VIEW_FLAGS = "-bhSu"
    def SAMTOOLS_SORT_FLAGS = SAMTOOLS_THREADS    

    transform(".fastq.gz") to(".bam") {
        exec """
            module load bwa &&
            module load samtools &&

            SAMPLE_NAME=\$(basename $output.prefix.prefix) &&

            PLATFORM="genomics" &&

            bwa mem $BWA_FLAGS -R \"@RG\\tID:\${SAMPLE_NAME}\\tSM:\${SAMPLE_NAME}\\tPL:illumina\\tLB:\${SAMPLE_NAME}\\tPU:\${PLATFORM}\" $BWA_REF $input | samtools view ${SAMTOOLS_VIEW_FLAGS} - | samtools sort ${SAMTOOLS_SORT_FLAGS} -T \${TMP}/\${SAMPLE_NAME} -  > ${output} &&

            samtools flagstat ${output} 1>&2
        ""","BWA_se"
    }
}

