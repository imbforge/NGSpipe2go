DE_DESeq2_MM = {
    doc title: "DE_DESeq2_MM",
        desc:  "Differential expression analysis using linear models and DESeq2, and considering multimappers",
        constraints: "Needs the results from dupRadar. Use the same GTF annotation as in subread_count",
        bpipe_version: "tested with bpipe 0.9.9.7",
        author: "Sergi Sayols"

    output.dir = DE_DESeq2_MM_vars.outdir
    def INPUT_READS_DIR = DE_DESeq2_MM_vars.cwd
    def DE_DESeq2_MM_FLAGS =
        (DE_DESeq2_MM_vars.targets   ? " targets="   + DE_DESeq2_MM_vars.targets   : "") +
        (DE_DESeq2_MM_vars.contrasts ? " contrasts=" + DE_DESeq2_MM_vars.contrasts : "") +
        (DE_DESeq2_MM_vars.filter    ? " filter="    + DE_DESeq2_MM_vars.filter    : "") +
        (DE_DESeq2_MM_vars.prefix    ? " prefix="    + DE_DESeq2_MM_vars.prefix    : "") +
        (DE_DESeq2_MM_vars.suffix    ? " suffix="    + DE_DESeq2_MM_vars.suffix    : "") +
        (DE_DESeq2_MM_vars.cwd       ? " cwd="       + DE_DESeq2_MM_vars.cwd       : "") +
        (DE_DESeq2_MM_vars.outdir    ? " out="       + DE_DESeq2_MM_vars.outdir    : "") +
        (DE_DESeq2_MM_vars.genes     ? " gtf="       + DE_DESeq2_MM_vars.genes     : "") +
        (DE_DESeq2_MM_vars.pattern   ? " pattern="   + DE_DESeq2_MM_vars.pattern      : "") +
        (DE_DESeq2_MM_vars.FDR       ? " FDR="       + DE_DESeq2_MM_vars.FDR        : "") +
        (DE_DESeq2_MM_vars.FC        ? " FC="        + DE_DESeq2_MM_vars.FC        : "") +
        (DE_DESeq2_MM_vars.extra     ? " "           + DE_DESeq2_MM_vars.extra     : "") 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // The DE_DESeq2_MM module is not using any of its inputs, but needs to check their
    // time stamp in order to know, if DE_DESeq2_MM should run (in case of pre-existing 
    // results). This can be done by outputting/echo'ing all inputs. In order to not 
    // confuse the pipeline user, this output is written to /dev/null
    // --- THE echo COMMAND BELOW MUST NOT BE REMOVED ---

    // run the chunk
    // should match deseq2.groovy, adding a step in between to convert all dupRadar input counts to HTSeq
    produce("DE_DESeq2.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            echo $inputs > /dev/null  &&

            if [[ ! -e "$INPUT_READS_DIR" ]]; then
                mkdir "$INPUT_READS_DIR";
            fi &&

            for f in $DE_DESeq2_MM_vars.dupradar_outdir/*.tsv; do
                F=\$(basename \$f) ;
                tail -n +2 $f | cut -f1,3 | sort -k1,1 > "$INPUT_READS_DIR/\${F%_dupRadar.tsv}.readcounts.tsv" ;
            done &&

            Rscript ${PIPELINE_ROOT}/tools/DE_DESeq2/DE_DESeq2.R $DE_DESeq2_MM_FLAGS
        ""","DE_DESeq2_MM"
    }
}

