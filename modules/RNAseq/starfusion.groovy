STAR_Fusion = {
    doc title: "STAR-Fusion",
        desc:  "detection of fusion transcripts from RNA-Seq data",
        constraints: "tab-delimited summary file identifying the fusion pairs. Works only with PE data",
        bpipe_version: "tested with bpipe 0.9.9",
        author: "Giuseppe Petrosino"

    output.dir = STAR_Fusion_vars.outdir

    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /.R1.fastq.gz/).replaceFirst("")

    def STARFUSION_FLAGS =
        (STAR_Fusion_vars.threads    ? " --CPU "            + STAR_Fusion_vars.threads    : "") + 
        (STAR_Fusion_vars.genome_lib ? " --genome_lib_dir " + STAR_Fusion_vars.genome_lib : "") +
        (STAR_Fusion_vars.extra      ? " "                  + STAR_Fusion_vars.extra      : "")

    def TOOL_ENV = prepare_tool_env("starfusion", tools["starfusion"]["version"], tools["starfusion"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(OUTPUTFILE + "_starfusion.done") {   // change it to whatever STAR-Fusion produces, and remove the touch $output, it's useless!
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            STAR-Fusion $STARFUSION_FLAGS --tmpdir \${TMP}/\$(basename $output.prefix) --left_fq $input1 --right_fq $input2 --output_dir $output.prefix &&
            touch $output
        ""","STAR_Fusion"
    }
}
