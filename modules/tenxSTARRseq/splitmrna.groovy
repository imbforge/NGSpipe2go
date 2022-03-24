SplitmRNA = {
        doc title: "SplitmRNA",
        desc:  "Split STARR-seq and endogenous mRNA 10X reads using cutadapt",
        constraints: "Requires paired end FASTQ files, with the STARR-seq mRNA-specific sequence at the 3' end of Read 1 (after the 10x barcode + UMI). Only supports compressed FASTQ files. Required naming scheme is bcl2fastq conventional '_S.*_L00?_R[12]_001.fastq.gz' pattern.",
        author: "Martin Oti, Antonio Domingues, Anke Busch, Nastasja Kreim, Frank Rühle"

    output.dir = SplitmRNA_vars.outdir

    def SAMPLENAME = new File(input.prefix.prefix)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_PRUNED = SAMPLENAME_BASE.replace(/_S.*_L00._R._001/, "") // delete Illumina fastq name suffix
    def SAMPLENAME_BASE_R2 = SAMPLENAME_BASE.replace("_R1_001", "_R2_001") // needed for paired end processing

    // create the log folder if it doesn't exists
    def SPLITMRNA_STATSDIR = new File( SplitmRNA_vars.statsdir)
    if (!SPLITMRNA_STATSDIR.exists()) {
        SPLITMRNA_STATSDIR.mkdirs()
    }

    def SPLITMRNA_FLAGS =
        (SplitmRNA_vars.adapter_sequence   ? " --adapter starr='"  + SplitmRNA_vars.adapter_sequence + "\$'"  : "") +
        (SplitmRNA_vars.minimum_overlap    ? " --overlap="         + SplitmRNA_vars.minimum_overlap           : "") +
        (Cutadapt_vars.minimum_length_keep ? " --minimum-length "  + Cutadapt_vars.minimum_length_keep        : "") +
        (Cutadapt_vars.maximum_length_keep ? " --maximum-length "  + Cutadapt_vars.maximum_length_keep        : "") +
        (SplitmRNA_vars.errorrate          ? " --error-rate "      + SplitmRNA_vars.errorrate                 : "") +
        (SplitmRNA_vars.action             ? " --action "          + SplitmRNA_vars.action                    : "") +
        (SplitmRNA_vars.extra              ? " "                   + SplitmRNA_vars.extra                     : "")

    def ADAPTERLENGTH = SplitmRNA_vars.adapter_sequence.length()
        
    def TOOL_ENV = prepare_tool_env("cutadapt", tools["cutadapt"]["version"], tools["cutadapt"]["runenv"])
    def PREAMBLE = get_preamble("cutadapt")

    transform(".fastq.gz") to ("_endogenous.fastq.gz", "_starr.fastq.gz") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            cutadapt $SPLITMRNA_FLAGS --output=\${TMP}/${SAMPLENAME_BASE}_{name}.fastq.gz --paired-output \${TMP}/${SAMPLENAME_BASE_R2}_{name}.fastq.gz --untrimmed-output=\${TMP}/${SAMPLENAME_BASE}_endogenous.fastq.gz --untrimmed-paired-output \${TMP}/${SAMPLENAME_BASE_R2}_endogenous.fastq.gz $input1 $input2 1> ${SPLITMRNA_STATSDIR}/${SAMPLENAME_BASE_PRUNED}.cutadapt_splitmrna.log &&

            mv -t $output.dir \${TMP}/${SAMPLENAME_BASE}*.fastq.gz    &&
            mv -t $output.dir \${TMP}/${SAMPLENAME_BASE_R2}*.fastq.gz
        ""","SplitmRNA"
    }
}


