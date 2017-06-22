DE_DESeq2 = {
    doc title: "DE_DESeq2",
        desc:  "Differential expression analysis using linear models and DESeq2",
        constraints: "Only simple contrasts in 1-factor design. Include always the intercept. Always gene filtering",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = DE_DESeq2_OUTDIR.replaceFirst("out=", "")
    def DE_DESeq2_FLAGS = DE_DESeq2_TARGETS   + " " + 
                          DE_DESeq2_CONTRASTS + " " +
                          DE_DESeq2_MMATRIX   + " " +
                          DE_DESeq2_FILTER    + " " +
                          DE_DESeq2_PREFIX    + " " +
                          DE_DESeq2_SUFFIX    + " " +
                          DE_DESeq2_CWD       + " " +
                          DE_DESeq2_OUTDIR    + " " +
                          DE_DESeq2_GENES     + " " +
                          DE_DESeq2_BASE      + " " +
                          DE_DESeq2_EXTRA

    // run the chunk
    produce("DE_DESeq2.RData") {
        exec """
            module load R/${R_VERSION} &&

            Rscript ${TOOL_DESeq2}/DE_DESeq2.R $DE_DESeq2_FLAGS
        ""","DE_DESeq2"
    }
}

