peak_annotation = {
    doc title: "peak_annotation",
        desc: "ChIP peak annotation, comparison and visualization using ChIPSeeker package",
        constraints: "",
        bpipe_version:"",
        author:"Sri Dewi"

    var subdir : ""
    output.dir = peak_annotation_vars.outdir + "/$subdir" 
 
    def PEAK_ANNOTATION_FLAGS =
        (peak_annotation_vars.files           ? " peakData="       + peak_annotation_vars.files + "/$subdir"  : "") +
        (peak_annotation_vars.transcript_type ? " transcriptType=" + peak_annotation_vars.transcript_type     : "") +
        (peak_annotation_vars.transcript_db   ? " transcriptDb="   + peak_annotation_vars.transcript_db       : "") +
        (peak_annotation_vars.orgdb           ? " orgDb="          + peak_annotation_vars.orgdb               : "") +
        (peak_annotation_vars.regiontss       ? " regionTSS="      + peak_annotation_vars.regiontss           : "") +
        (peak_annotation_vars.targets         ? " targets="        + peak_annotation_vars.targets             : "") +
        (peak_annotation_vars.orderby         ? " orderby="        + peak_annotation_vars.orderby             : "") +
        (peak_annotation_vars.outdir          ? " out="            + peak_annotation_vars.outdir + "/$subdir" : "") +
        (peak_annotation_vars.extra           ? " "                + peak_annotation_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce("Peak_Annotation.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/Peak_Annotation/Peak_Annotation.R $PEAK_ANNOTATION_FLAGS;
        ""","peak_annotation"
    }
}
