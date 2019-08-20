// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/starfusion.vars.groovy"

STAR_Fusion = {
    doc title: "STAR-Fusion",
        desc:  "detection of fusion transcripts from RNA-Seq data",
        constraints: "tab-delimited summary file identifying the fusion pairs. Works only with PE data",
        bpipe_version: "tested with bpipe 0.9.9",
        author: "Giuseppe Petrosino"

    output.dir = STARFUSION_OUTDIR

    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")

    def STARFUSION_FLAGS = STARFUSION_THREADS + " " +
                           STARFUSION_GENOME_LIB

    def TOOL_ENV = prepare_tool_env("starfusion", tools["starfusion"]["version"], tools["starfusion"]["runenv"])
    def PREAMBLE = get_preamble("STAR_Fusion")

    produce(OUTPUTFILE + "_starfusion.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            STAR-Fusion $STARFUSION_FLAGS --tmpdir \${TMP}/\$(basename $output.prefix) --left_fq $input1 --right_fq $input2 --output_dir $output.prefix;
            touch $output
        ""","STAR_Fusion"
    }
}
