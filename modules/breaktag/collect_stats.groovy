collect_stats = {
    doc title: "collect stats",
    desc:  "collect breaktag DSB stats",
    constraints: "none",
    author: "Sergi Sayols"

    output.dir = collect_stats_vars.outdir

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".strandless.bed.gz") to(".txt") {
        exec """
            ${PREAMBLE} &&
            zcat $input | awk '{i+=\$5} END {print "breaks: ", i; print "loci: ", NR;}' > $output
        """
    }
    forward input
}

