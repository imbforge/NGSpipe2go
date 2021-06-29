BWA_pe = {
    doc title: "BWA PE alignment",
    desc:  "Align paired end reads",
    constraints: "Only works with compressed input. Set all global vars. Samtools version >= 1.2.",
    bpipe_version: "tested with bpipe 0.9.8.7",
    author: "Oliver Drechsel, Anke Busch"

    output.dir = BWA_vars.outdir

    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    def BWA_FLAGS = "-M " +
        (BWA_vars.threads ? "-t " + BWA_vars.threads : "" ) +
        (BWA_vars.extra   ?         BWA_vars.extra   : "" )

    def SAMTOOLS_VIEW_FLAGS = "-bhSu"
    def SAMTOOLS_SORT_FLAGS = 
        (BWA_vars.samtools_threads ? " -@ " + BWA_vars.samtools_threads : "" )

    def TOOL_ENV = prepare_tool_env("bwa", tools["bwa"]["version"], tools["bwa"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(OUTPUTFILE + ".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bwa mem $BWA_FLAGS -R \"@RG\\tID:${OUTPUTFILE}\\tSM:${OUTPUTFILE}\\tPL:illumina\\tLB:${OUTPUTFILE}\\tPU:genomics\" $BWA_vars.ref $input1 $input2 | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/${OUTPUTFILE}_sort -  > ${output} &&

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

    output.dir = BWA_vars.outdir

    File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    def BWA_FLAGS = "-M " +
        (BWA_vars.threads ? "-t " + BWA_vars.threads : "" ) +
        (BWA_vars.extra   ?         BWA_vars.extra   : "" )

    def SAMTOOLS_VIEW_FLAGS = "-bhSu"
    def SAMTOOLS_SORT_FLAGS = 
        (BWA_vars.samtools_threads ? " -@ " + BWA_vars.samtools_threads : "" )

    def TOOL_ENV = prepare_tool_env("bwa", tools["bwa"]["version"], tools["bwa"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".fastq.gz") to(".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bwa mem $BWA_FLAGS -R \"@RG\\tID:\${OUTPUTFILE}\\tSM:\${OUTPUTFILE}\\tPL:illumina\\tLB:\${OUTPUTFILE}\\tPU:genomics\" $BWA_vars.ref $input | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/${OUTPUTFILE}_sort - > ${output} &&

            samtools flagstat ${output} 1>&2
        ""","BWA_se"
    }
}

