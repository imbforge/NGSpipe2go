load MODULE_FOLDER + "SmallRNAseq/aggregate_mapped_counts.vars.groovy"

AggregateMappedCounts = {
    doc title: "AggregateMappedCounts",
    desc: "Puts reads counts from all the bams in a single file read to be read in R as normalization factors.",
    constraints: "none",
    author: "António Domingues"

    output.dir = FINAL_COUNTS_OUT

    from("*.mappedreads.txt") produce("mappedReads.txt"){

        exec """
            cat $inputs > $output
      ""","AggregateMappedCounts"
    }
}
