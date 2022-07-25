FastqScreen = {
    doc title: "FastqScreen",
        desc:  "Quality control of input file against various contaminants",
        constraints: "Only supports compressed FASTQ files",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim, modified by Frank Ruehle"

    var subdir : ""
    output.dir = FastqScreen_vars.outdir + "/$subdir"
    def FASTQSCREEN_FLAGS = 
        (FastqScreen_vars.threads ? " --threads " + FastqScreen_vars.threads : "") +
        (FastqScreen_vars.extra   ? " "           + FastqScreen_vars.extra   : "")

    def TOOL_ENV = prepare_tool_env("fastqscreen", tools["fastqscreen"]["version"], tools["fastqscreen"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform("*.fastq.gz") to("_fastqscreen.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [ ! -e "$output.prefix" ]; then
                mkdir $output.prefix;
            fi &&
            fastqreference=${FastqScreen_vars.genome}${FastqScreen_vars.conf};
            references=(\${fastqreference//,/ });
            for i in "\${!references[@]}"; do
                reference=(\${references[i]//::/ });
                echo -e "DATABASE\t\${reference[0]}\t\${reference[1]}" >> $output.prefix/fastqscreen.conf;
            done;
            fastq_screen $FASTQSCREEN_FLAGS --conf $output.prefix/fastqscreen.conf --outdir $output.prefix $inputs;
            touch $outputs
        ""","FastqScreen"
    }

    forward inputs
}

