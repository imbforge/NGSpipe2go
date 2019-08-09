// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/peak_annotation.vars.groovy"

peak_annotation = {

    doc title: "peak_annotation",
        desc: "ChIP peak annotation, comparison and visualization using ChIPSeeker package",
        constraints: "",
        bpipe_version:"",
        author:"Sri Dewi"

    output.dir=peak_annotation_OUTDIR.replaceFirst("out=","")

    def peak_annotation_FLAGS = peak_annotation_FILES + " " +
        peak_annotation_TRANSCRIPT_TYPE + " " + 
        peak_annotation_TRANSCRIPT_DB + " " + 
        peak_annotation_ORGDB + " " + 
        peak_annotation_REGIONTSS + " " + 
        peak_annotation_OUTDIR + " " +
        peak_annotation_EXTRA

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    produce("Peak_Annotation.RData") {
        exec """
            ${TOOL_ENV} &&

            Rscript ${PIPELINE_ROOT}/tools/Peak_Annotation/Peak_Annotation.R $peak_annotation_FLAGS;
        ""","peak_annotation"
    }
}
