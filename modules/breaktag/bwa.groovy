bwa = {
    doc title: "BWA SR/PE alignment",
    desc:  "Align SR/PE reads with BWA",
    constraints: "none",
    author: "Sergi Sayols"

    output.dir = BWA_vars.outdir

    def File f = new File(input1)
    def OUTFILE = (f.getName() =~ /(.R1)*.filtered.fastq.gz/).replaceFirst(".bam")

    def BWA_INPUT = (BWA_vars.paired ? "$input1 $input2" : "$input")
    def BWA_FLAGS =
        (BWA_vars.threads ? " -t " + BWA_vars.threads : "" ) +
        (BWA_vars.extra   ? " "    + BWA_vars.extra   : "" )

    def SAMTOOLS_VIEW_FLAGS = "-bhSu" +
        (BWA_vars.minqual ? " -q " + BWA_vars.minqual : "")
    def SAMTOOLS_SORT_FLAGS = 
        (BWA_vars.samtools_threads ? " -@ " + BWA_vars.samtools_threads : "" )

    def TOOL_ENV = prepare_tool_env("bwa", tools["bwa"]["version"], tools["bwa"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(OUTFILE) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bwa mem $BWA_FLAGS $BWA_vars.ref $BWA_INPUT | \
              samtools view $SAMTOOLS_VIEW_FLAGS - | \
              samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/${OUTPUTFILE}_sort - > ${output} &&

            samtools index ${output}
            ""","BWA_pe"
    }
}

