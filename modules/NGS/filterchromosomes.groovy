FilterChr = {
    doc title: "FilterChr",
        desc:  "When mapping to full genome, including unassembled contigs, remove those extra contiguous before proceeding for further analysis. The goal is to increase speed and decrease disk space usage. Source: https://www.biostars.org/p/171791/#171819",
        constraints: "Requires a file with the list of chromosomes to keep.",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "António Domingues"

    output.dir=FilterChr_vars.outdir

    def PREAMBLE = get_preamble("FilterChr")

    transform(".bam") to (".chrOnly.bam") {
        exec """
            ${PREAMBLE} &&

            chroms=\$(cut -f1 ${FilterChr_vars.file}) && 
            samtools view -@ ${FilterChr_vars.threads} -b $input \$chroms > $output &&
            samtools index $output
        ""","FilterChr"
    }
}

