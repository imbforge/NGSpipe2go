load MODULE_FOLDER + "ChIPseq/peak_annotation.vars.groovy"

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


    produce("Peak_Annotation.RData") {
        exec """
            module load R/${R_VERSION} &&

            Rscript ${TOOL_PEAK_ANNOTATION}/Peak_Annotation.R $peak_annotation_FLAGS;
        ""","peak_annotation"
    }
}
