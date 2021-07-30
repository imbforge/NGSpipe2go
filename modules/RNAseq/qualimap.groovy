qualimap = {
    doc title: "Qualimap",
        desc:  "Call qualimap to do rnaseq qualitycontrol",
        author: "Nastasja Kreim"

    output.dir = qualimap_vars.outdir
    // no|yes|reverse
    if(qualimap_vars.stranded == "no") {
        qualimap_vars.protocol = "non-strand-specific"
    }
    else if (qualimap_vars.stranded == "yes") {
        qualimap_vars.protocol = "strand-specific-forward"
    }
    else {
        qualimap_vars.protocol = "strand-specific-reverse"
    }
    if(qualimap_vars.paired){
        qualimap_vars.extra = qualimap_vars.extra + " -pe"
    }

    def TOOL_ENV = prepare_tool_env("qualimap", tools["qualimap"]["version"], tools["qualimap"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bam") to("_counts.txt") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            unset DISPLAY;
            echo $output.prefix;
            qualimap rnaseq -bam $input -outdir ${output.prefix}_qualimap -outformat html -gtf ${qualimap_vars.genesgtf} -oc $output -p ${qualimap_vars.protocol} ${qualimap_vars.extra}
        ""","qualimap"
    }

    forward input
}

