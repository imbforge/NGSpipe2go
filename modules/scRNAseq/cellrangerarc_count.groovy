cellrangerarc_count = {
    doc title: "Cell Ranger ARC (10x Multiome) alignment",
        desc:  "Align single/paired end reads",
        constraints: "Paired end reads expected to have _S1_L001_R[12]_001.fastq.gz suffix.",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle, Martin Oti"

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_PRUNED = SAMPLENAME_BASE.replaceAll(~/_gex_S[\d]+_L[\d]+_.*/, "").replaceAll(~/_atac_S[\d]+_L[\d]+_.*/, "")
    def SAMPLENAME_GEX = SAMPLENAME_BASE.replaceAll("_atac_", "_gex_")
    def SAMPLENAME_ATAC = SAMPLENAME_BASE.replaceAll("_gex_", "_atac_")
    def SAMPLEDIR = cellrangerarc_count_vars.fastqdir

    output.dir = cellrangerarc_count_vars.outdir + "/" 

    // cellranger flags
    def CELLRANGERARC_FLAGS =
        " --id="     + SAMPLENAME_BASE +
        " --reference=" + cellrangerarc_count_vars.reference + 
        (cellrangerarc_count_vars.cores        ? " --localcores="    + cellrangerarc_count_vars.cores        : "") + 
        (cellrangerarc_count_vars.localmem     ? " --localmem="      + cellrangerarc_count_vars.localmem     : "") + 
        (cellrangerarc_count_vars.extra        ? " "                 + cellrangerarc_count_vars.extra        : "")

    def TOOL_ENV = prepare_tool_env("cellrangerarc", tools["cellrangerarc"]["version"], tools["cellrangerarc"]["runenv"]) 
    def PREAMBLE = get_preamble("cellrangerarc")

    transform('.fastq.gz') to('.bam') {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            echo "fastqs,sample,library_type" > libraries_${SAMPLENAME_BASE}.csv &&
            echo "${SAMPLEDIR},${SAMPLENAME_PRUNED}_gex,Gene Expression" >> libraries_${SAMPLENAME_BASE}.csv &&
            echo "${SAMPLEDIR},${SAMPLENAME_PRUNED}_atac,Chromatin Accessibility" >> libraries_${SAMPLENAME_BASE}.csv &&

            cellranger-arc count --libraries=libraries_${SAMPLENAME_BASE}.csv $CELLRANGERARC_FLAGS &&
            mv ${SAMPLENAME_BASE}/outs/gex_possorted_bam.bam ${output.dir}/${SAMPLENAME_GEX}.bam &&
            mv ${SAMPLENAME_BASE}/outs/gex_possorted_bam.bam.bai ${output.dir}/${SAMPLENAME_GEX}.bam.bai &&
            mv ${SAMPLENAME_BASE}/outs/atac_possorted_bam.bam ${output.dir}/${SAMPLENAME_ATAC}.bam &&
            mv ${SAMPLENAME_BASE}/outs/atac_possorted_bam.bam.bai ${output.dir}/${SAMPLENAME_ATAC}.bam.bai &&
            mv libraries_${SAMPLENAME_BASE}.csv $output.dir &&
            mv ${SAMPLENAME_BASE} $output.dir

        ""","cellrangerarc_count"
    }
}


