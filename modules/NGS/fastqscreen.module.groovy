// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/fastqscreen.vars.groovy"

FastqScreen = {
    doc title: "FastScreen",
        desc:  "Quality control of input file against various contaminants",
        constraints: "Only supports compressed FASTQ files",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = FASTQSCREEN_OUTDIR
    def FASTQSCREEN_FLAGS = "-threads" + FASTQSCREEN_THREADS + " " + FASTQSCREEN_PARAM

    def TOOL_ENV = prepare_tool_env("fastqscreen", tools["fastqscreen"]["version"], tools["fastqscreen"]["runenv"])

    transform(".fastq.gz") to("_fastqscreen.done") {
        exec """
            ${TOOL_ENV} &&

            if [ ! -e "$output.prefix" ]; then
                mkdir $output.prefix;
            fi &&
            fastqreference=$FASTQSCREEN_CONF;
            references=(\${fastqreference//,/ });
            for i in "\${!references[@]}"; do
                reference=(\${references[i]//::/ });
                echo -e "DATABASE\t\${reference[0]}\t\${reference[1]}" >> $output.prefix/fastqscreen.conf;
            done;
            fastq_screen $FASTQSCREEN_PARAM --conf $output.prefix/fastqscreen.conf $FASTQSCREEN_PARAM --outdir $output.prefix $input;
            touch $output
        ""","FastqScreen"
    }

    forward input
}

