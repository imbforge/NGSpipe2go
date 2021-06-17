cellranger_count = {
    doc title: "Cellranger alignment",
        desc:  "Align single/paired end reads",
        constraints: "Paired end reads expected to have _S1_L001_R[12]_001.fastq.gz suffix.",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLEDIR = SAMPLENAME.getParentFile()
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_PRUNED = SAMPLENAME_BASE.replaceAll(~/_S[\d]+_L[\d]+_.*/, "") 

    output.dir = cellranger_count_vars.outdir + "/" 

    // cellranger flags
    def CELLRANGER_FLAGS =
        " --id="     + SAMPLENAME_BASE +
        " --fastqs=" + SAMPLEDIR + 
        " --sample=" + SAMPLENAME_PRUNED + 
        " --transcriptome=" + cellranger_count_vars.transcriptome + 
        (cellranger_count_vars.expect_cells ? " --expect-cells="  + cellranger_count_vars.expect_cells : "") + 
        (cellranger_count_vars.cores        ? " --localcores="    + cellranger_count_vars.cores        : "") + 
        (cellranger_count_vars.localmem     ? " --localmem="      + cellranger_count_vars.localmem     : "") + 
        (cellranger_count_vars.extra        ? " "                 + cellranger_count_vars.extra        : "")

    def TOOL_ENV = prepare_tool_env("cellranger", tools["cellranger"]["version"], tools["cellranger"]["runenv"]) 
    def PREAMBLE = get_preamble("cellranger")

    transform('.fastq.gz') to('.bam') {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            cellranger count $CELLRANGER_FLAGS &&
            mv ${SAMPLENAME_BASE}/outs/possorted_genome_bam.bam $output && 
            mv ${SAMPLENAME_BASE}/outs/possorted_genome_bam.bam.bai ${output.dir}/${SAMPLENAME_BASE}.bam.bai &&
            mv ${SAMPLENAME_BASE} $output.dir 

        ""","cellranger_count"
    }
}


