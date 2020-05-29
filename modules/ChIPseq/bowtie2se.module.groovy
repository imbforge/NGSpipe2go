// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/bowtie2se.vars.groovy"

bowtie2_se = {
    doc title: "Bowtie2 SE alignment",
        desc:  "Align single reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = bowtie2_se_vars.mapped

    def BOWTIE2_FLAGS = "-q " +
        (bowtie2_se_vars.quals   ?  " "     + bowtie2_se_vars.quals   : "") + 
        (bowtie2_se_vars.mm_seed ?  " "     + bowtie2_se_vars.mm_seed : "") + 
        (bowtie2_se_vars.insert  ?  " -L "  + bowtie2_se_vars.insert  : "") +
        (bowtie2_se_vars.threads ?  " -p "  + bowtie2_se_vars.threads : "") + 
        (bowtie2_se_vars.extra   ?  " "     + bowtie2_se_vars.extra   : "") +
        (bowtie2_se_vars.ref     ?  " -x "  + bowtie2_se_vars.ref     : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie2_se_vars.samtools_threads ? " -@ " + bowtie2_se_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie2", tools["bowtie2"]["version"], tools["bowtie2"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie2_se")

    transform(".fastq.gz") to (".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            bowtie2 $BOWTIE2_FLAGS -U $input | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output;
        ""","bowtie2_se"
    }
}



