// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/qualimap.vars.groovy"

qualimap = {
    doc title: "Qualimap",
        desc:  "Call qualimap to do rnaseq qualitycontrol",
        author: "Nastasja Kreim"

    output.dir = QUALIMAP_OUTDIR
    // no|yes|reverse
    if(QUALIMAP_STRANDED == "no") {
        QUALIMAP_PROTOCOL = "non-strand-specific"
    }
    else if (QUALIMAP_STRANDED == "yes") {
        QUALIMAP_PROTOCOL = "strand-specific-forward"
    }
    else {
        QUALIMAP_PROTOCOL = "strand-specific-reverse"
    }
    if(QUALIMAP_PAIRED == "yes"){
        QUALIMAP_EXTRA = QUALIMAP_EXTRA + " -pe"
    }

    def TOOL_ENV = prepare_tool_env("qualimap", tools["qualimap"]["version"], tools["qualimap"]["runenv"])
    def PREAMBLE = get_preamble("qualimap")

    transform(".bam") to("_counts.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            unset DISPLAY;
            echo $output.prefix;
            qualimap rnaseq -bam $input -outdir ${output.prefix}_qualimap -outformat html $QUALIMAP_GENESGTF -oc $output -p $QUALIMAP_PROTOCOL $QUALIMAP_EXTRA
        ""","qualimap"
    }

    forward input
}

