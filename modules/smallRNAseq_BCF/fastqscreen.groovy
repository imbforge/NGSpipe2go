FastQScreen = {
    doc title: "FastQScreen",
    desc:  "Quality control of input file against various contaminants",
    constraints: "Only supports compressed FASTQ files",
    author: "Nastasja Kreim, Anke Busch"

    output.dir = FastQScreen_vars.outdir
    def FASTQSCREEN_FLAGS =
        (FastQScreen_vars.threads ? " --threads " + FastQScreen_vars.threads : "") +
        (FastQScreen_vars.param   ? " "           + FastQScreen_vars.param   : "")

    def TOOL_ENV = prepare_tool_env("fastqscreen", tools["fastqscreen"]["version"], tools["fastqscreen"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:"", input:new File(input1.prefix).getName())

    transform(".fastq.gz") to ("_fastqscreen.done") {
        def SAMPLENAME = output.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [ ! -e "$output.prefix" ]; then
                mkdir $output.prefix;
            fi &&

            FASTQREFERENCES=$FastQScreen_vars.conf;
            REFERENCES=(\${FASTQREFERENCES//,/ });

            for i in "\${!REFERENCES[@]}"; do
                REFERENCE=(\${REFERENCES[i]//::/ });
                echo -e "DATABASE\t\${REFERENCE[0]}\t\${REFERENCE[1]}" >> $output.prefix/fastqscreen.conf;
            done &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            fastq_screen $FASTQSCREEN_FLAGS --conf $output.prefix/fastqscreen.conf --outdir $output.prefix $input 2> $output.dir/\${SAMPLENAME_BASE}.log &&
            touch $output
        ""","FastQScreen"
    }
    forward input
}

