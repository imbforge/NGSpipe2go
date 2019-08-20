// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/blacklist_filter.vars.groovy"

blacklist_filter = {
    doc title: "blacklist_filter",
        desc: "Remove peaks overlapping blacklisted genomic regions",
        constraints: "",
        bpipe_version:"",
        author:"Giuseppe Petrosino"

    output.dir=blacklist_filter_OUTDIR.replaceFirst("out=", "")
    def blacklist_filter_FLAGS = blacklist_filter_FILES    + " " +
                                 blacklist_filter_BED_FILE + " " + 
                                 blacklist_filter_OUTDIR   + " " +
                                 blacklist_filter_EXTRA

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble("blacklist_filter")

    produce("BlackList_Filter.RData") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            Rscript ${PIPELINE_ROOT}/tools/BlackList_Filter/BlackList_Filter.R $blacklist_filter_FLAGS;
        ""","blacklist_filter"
    }
}
