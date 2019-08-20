// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/fastqscreen.vars.groovy"

FastQScreen = {
    doc title: "FastQScreen",
    desc:  "Quality control of input file against various contaminants",
    constraints: "Only supports compressed FASTQ files",
    author: "Nastasja Kreim, Anke Busch"

    output.dir   = FASTQSCREEN_OUTDIR
    def FASTQSCREEN_FLAGS = "--threads " + FASTQSCREEN_THREADS + " " + FASTQSCREEN_PARAM

    def TOOL_ENV = prepare_tool_env("fastq_screen", tools["fastq_screen"]["version"], tools["fastq_screen"]["runenv"])
    def PREAMBLE = get_preamble("FastQScreen")

    transform(".fastq.gz") to ("_fastqscreen.done") {
        def SAMPLENAME = output.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [ ! -e "$output.prefix" ]; then
                mkdir $output.prefix;
            fi &&

            FASTQREFERENCES=$FASTQSCREEN_CONF;
            REFERENCES=(\${FASTQREFERENCES//,/ });

            for i in "\${!REFERENCES[@]}"; do
                REFERENCE=(\${REFERENCES[i]//::/ });
                echo -e "DATABASE\t\${REFERENCE[0]}\t\${REFERENCE[1]}" >> $output.prefix/fastqscreen.conf;
            done;

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            fastq_screen $FASTQSCREEN_FLAGS --conf $output.prefix/fastqscreen.conf --outdir $output.prefix $input 2> $output.dir/\${SAMPLENAME_BASE}.log;
            touch $output
        ""","FastQScreen"
    }
    forward input
}

