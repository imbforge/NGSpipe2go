collect_stats = {
    doc title: "collect stats",
    desc:  "collect breaktag DSB stats",
    constraints: "none",
    author: "Sergi Sayols"

    output.dir = collect_stats_vars.outdir

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform("bed") to("txt") {
        exec """
            ${PREAMBLE} &&
            echo "$input stats" > $output
        """
    }
    forward input
}

