filbowtie2unique = {
    doc title: "filter out multimapping reads from bowtie2 out and removes duplicate reads, if the parameter is set",
        desc:  "filter out multimapping reads from bowtie2 output bam file and removes PCR duplicates if the parameter is set",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Nastasja Kreim, Sivarajan Karunanithi"

    output.dir = filbowtie2unique_vars.mapped

    def FILBOWTIE2_FLAGS = (filbowtie2unique_vars.paired ? " -f 2 -q $filbowtie2unique_vars.samtools_mapq_pe" : " -F 4 -q $filbowtie2unique_vars.samtools_mapq_se")

    def TOOL_ENV = prepare_tool_env("bamutil", tools["bamutil"]["version"], tools["bamutil"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())
    // This variable allows us to create flexible output file names depending on the parameters used
    def OUTBAMFILE = (dupremoval_vars.remove_pcr_dups ? ".unique.duprm.bam" : ".unique.bam")

    transform(".bam") to (OUTBAMFILE) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            samtools view $FILBOWTIE2_FLAGS -bhu ${input} | samtools sort -@ $filbowtie2unique_vars.samtools_threads -T \${TMP}/\$(basename $output.prefix) -o ${input.prefix}_mmremoved.bam - &&
            if [[ "${dupremoval_vars.remove_pcr_dups}" == "true" ]]; then
                echo "Removing PCR duplicates";
                bam dedup --in ${input.prefix}_mmremoved.bam --out ${output} --log ${LOGS}/filbowtie2unique/\$(basename ${output.prefix})_dupmetrics.log --noPhoneHome --force --rmDups;
                rm ${input.prefix}_mmremoved.bam;
            else
                echo "Not removing the PCR duplicates";
                mv ${input.prefix}_mmremoved.bam ${output};
            fi;

        ""","filbowtie2unique"
    }
}

