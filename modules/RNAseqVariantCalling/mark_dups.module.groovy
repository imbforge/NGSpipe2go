// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/mark_dups.vars.groovy"

MarkDups = {
    doc title: "MarkDups",
        desc: "Call picard tools to mark with/without removing duplicated reads from a bam file",
        constraints: "Picard tools version >= 1.141"
        author: "Sergi Sayols, modified by Antonio Domingues"

    output.dir = OUTDIR_STAR2ND

    def JAVA_FLAGS = "-Xmx" + MARKDUPS_MAXMEM
    def MARKDUPS_FLAGS  = " REMOVE_DUPLICATES=" + MARKDUPS_REMOVE +
                          " CREATE_INDEX=" + MARKDUPS_INDEX +
                          " VALIDATION_STRINGENCY=" + MARKDUPS_VALIDATION +
                          " ASSUME_SORTED=TRUE"

    def TOOL_ENV = prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])

    transform(".rg.bam") to (".rg.duprm.bam"){
        exec """
            ${TOOL_ENV} &&

            echo 'VERSION INFO'  1>&2 &&
            echo \$(java -jar \${picard} MarkDuplicates --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            java $JAVA_FLAGS -jar \${picard} MarkDuplicates $MARKDUPS_FLAGS I=$input O=$output M=${input.prefix}_dupmetrics.tsv
        ""","MarkDups"
    }
}
