// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/deseq2_mm.vars.groovy"

DE_DESeq2_MM = {
    doc title: "DE_DESeq2_MM",
        desc:  "Differential expression analysis using linear models and DESeq2, and considering multimappers",
        constraints: "Needs the results from dupRadar. Use the same GTF annotation as in subread_count",
        bpipe_version: "tested with bpipe 0.9.9.7",
        author: "Sergi Sayols"

    output.dir = DE_DESeq2_MM_OUTDIR.replaceFirst("out=", "")
    def INPUT_READS_DIR = DE_DESeq2_MM_CWD.replaceFirst("cwd=", "")
    def DE_DESeq2_MM_FLAGS = DE_DESeq2_MM_TARGETS   + " " + 
                             DE_DESeq2_MM_CONTRASTS + " " +
                             DE_DESeq2_MM_MMATRIX   + " " +
                             DE_DESeq2_MM_FILTER    + " " +
                             DE_DESeq2_MM_PREFIX    + " " +
                             DE_DESeq2_MM_SUFFIX    + " " +
                             DE_DESeq2_MM_CWD       + " " +
                             DE_DESeq2_MM_OUTDIR    + " " +
                             DE_DESeq2_MM_GENES     + " " +
                             DE_DESeq2_MM_EXTRA

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("DE_DESeq2_MM")

    // run the chunk
    // should match deseq2.module.groovy, adding a step in between to convert all dupRadar input counts to HTSeq
    produce("DE_DESeq2.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [[ ! -e "$INPUT_READS_DIR" ]]; then
                mkdir "$INPUT_READS_DIR";
            fi &&

            for f in $DE_DESeq2_MM_DUPRADAR_OUTDIR/*.tsv; do
                F=\$(basename \$f) ;
                tail -n +2 $f | cut -f1,3 | sort -k1,1 > "$INPUT_READS_DIR/\${F%_dupRadar.tsv}.readcounts.tsv" ;
            done &&

            Rscript ${PIPELINE_ROOT}/tools/DE_DESeq2/DE_DESeq2.R $DE_DESeq2_MM_FLAGS
        ""","DE_DESeq2_MM"
    }
}

