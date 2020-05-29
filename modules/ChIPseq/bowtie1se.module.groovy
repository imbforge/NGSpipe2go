// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/bowtie1se.vars.groovy"

bowtie_se = {
    doc title: "Bowtie SE alignment",
        desc:  "Align single end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=1.2).",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = bowtie_se_vars.mapped

    def BOWTIE_FLAGS = "-q --sam "  +
        (bowtie_se_vars.quals   ? " "    + bowtie_se_vars.quals    : "") +
        (bowtie_se_vars.best    ? " --best --strata --tryhard --chunkmbs 256" : "") +
        (bowtie_se_vars.mm_seed ? " -n " + bowtie_se_vars.mm_seed  : "") +
        (bowtie_se_vars.insert  ? " -l " + bowtie_se_vars.insert   : "") +
        (bowtie_se_vars.maqerr  ? " -e " + bowtie_se_vars.maqerr   : "") +
        (bowtie_se_vars.multimap_mode && bowtie_se_vars.multimap ?
          (bowtie_se_vars.multimap_mode == "discard" ? " -m " : " -M ") + bowtie_se_vars.multimap : "") +
        (bowtie_se_vars.threads ? " -p " + bowtie_se_vars.threads  : "") +
        (bowtie_se_vars.extra   ? " "    + bowtie_se_vars.extra    : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie_se_vars.samtools_threads ? " -@ " + bowtie_se_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("bowtie_se")

    transform(".fastq.gz") to (".bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            zcat $input | bowtie $BOWTIE_FLAGS ${bowtie_se_vars.ref} - | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix)_bowtie1_sort - > $output
        ""","bowtie_se"
    }
}

