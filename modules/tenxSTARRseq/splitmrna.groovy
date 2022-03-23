SplitmRNA = {
        doc title: "SplitmRNA",
        desc:  "Split STARR-seq and endogenous mRNA 10X reads using cutadapt",
        constraints: "Requires paired end FASTQ files, with the STARR-seq mRNA-specific sequence at the 3' end of Read 1 (after the 10x barcode + UMI). Only supports compressed FASTQ files. Required naming scheme is bcl2fastq conventional '_S.*_L00?_R[12]_001.fastq.gz' pattern.",
        author: "Martin Oti, Antonio Domingues, Anke Busch, Nastasja Kreim, Frank Rühle"

    output.dir = SplitmRNA_vars.outdir

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_PRUNED = SAMPLENAME_BASE.replace("_S.*_L00._R._001", "") // delete .R1 in combined log file of pe design
    def SAMPLENAME_BASE_R2 = SAMPLENAME_BASE.replace("_R1_001", "_R2_001") // needed for paired end processing

    // create the log folder if it doesn't exists
    def SPLITMRNA_STATSDIR = new File( SplitmRNA_vars.statsdir)
    if (!SPLITMRNA_STATSDIR.exists()) {
        SPLITMRNA_STATSDIR.mkdirs()
    }

    def SPLITMRNA_FLAGS =
        (SplitmRNA_vars.adapter_sequence  ? " --adapter starr='"  + SplitmRNA_vars.adapter_sequence + "\$'"  : "") +
        (SplitmRNA_vars.minimum_overlap   ? " --overlap="         + SplitmRNA_vars.minimum_overlap           : "") +
        (SplitmRNA_vars.errorrate         ? " --error-rate "      + SplitmRNA_vars.errorrate                 : "") +
        (SplitmRNA_vars.extra             ? " "                   + SplitmRNA_vars.extra                     : "")

    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform("*.fastq.gz") to (".endogenous.fastq.gz", ".starr.fastq.gz") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            cutadapt $SPLITMRNA_FLAGS --output=\${TMP}/${SAMPLENAME_BASE}.{name}.fastq.gz --paired-output \${TMP}/${SAMPLENAME_BASE_R2}.{name}.fastq.gz --untrimmed-output=\${TMP}/${SAMPLENAME_BASE}.endogenous.fastq.gz --untrimmed-paired-output \${TMP}/${SAMPLENAME_BASE_R2}.endogenous.fastq.gz $input1 $input2.optional 1> ${SPLITMRNA_STATSDIR}/${SAMPLENAME_BASE_PRUNED}.cutadapt_splitmrna.log &&
		
            mv -t $output.dir \${TMP}/${SAMPLENAME_BASE_PRUNED}*.fastq.gz
        ""","SplitmRNA"
    }
}


