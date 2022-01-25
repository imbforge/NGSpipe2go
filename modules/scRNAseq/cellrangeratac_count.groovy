cellrangeratac_count = {
    doc title: "Cell Ranger ATAC alignment",
        desc:  "Align paired end reads with 10X barcode in i5 index",
        constraints: "Read 1, i5 & read 2 files expected to have _S1_L001_R[123]_001.fastq.gz suffix respectively. For this, run bcl2fastq with extra argument '--use-bases-mask Y*,I8,Y*,Y*'.",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle, Martin Oti"

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLEDIR = SAMPLENAME.getParentFile()
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_PRUNED = SAMPLENAME_BASE.replaceAll(~/_S[\d]+_L[\d]+_.*/, "") 

    output.dir = cellrangeratac_count_vars.outdir + "/" 

    // cellranger-atac count flags
    def CELLRANGERATAC_FLAGS =
        " --id="     + SAMPLENAME_BASE +
        " --fastqs=" + SAMPLEDIR + 
        " --sample=" + SAMPLENAME_PRUNED + 
        " --reference=" + cellrangeratac_count_vars.reference + 
        (cellrangeratac_count_vars.cores     ? " --localcores="  + cellrangeratac_count_vars.cores        : "") + 
        (cellrangeratac_count_vars.localmem  ? " --localmem="    + cellrangeratac_count_vars.localmem     : "") + 
        (cellrangeratac_count_vars.extra     ? " "               + cellrangeratac_count_vars.extra        : "")

    def TOOL_ENV = prepare_tool_env("cellrangeratac", tools["cellrangeratac"]["version"], tools["cellrangeratac"]["runenv"]) 
    def PREAMBLE = get_preamble("cellrangeratac")

    transform('.fastq.gz') to('.bam') {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            cellranger-atac count $CELLRANGERATAC_FLAGS &&
            mv ${SAMPLENAME_BASE}/outs/possorted_bam.bam $output && 
            mv ${SAMPLENAME_BASE}/outs/possorted_bam.bam.bai ${output.dir}/${SAMPLENAME_BASE}.bam.bai &&
            mv ${SAMPLENAME_BASE} $output.dir 

        ""","cellrangeratac_count"
    }
}


