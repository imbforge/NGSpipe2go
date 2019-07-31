load MODULE_FOLDER + "NGS/markdups2.vars.groovy"

MarkDups2 = {
    doc title: "MarkDups2",
        desc:  "Call bamUtil dedup tool to mark with/without removing duplicated reads from a bam file",
        constraints: "bamUtil tool version >= 1.0.13",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Giuseppe Petrosino"

    output.dir=MAPPED

    transform(".bam") to (".dupmarked.bam") {
        exec """
            module load bamUtil/${BAMUTIL_VERSION} &&

            bam dedup --in $input --out $output --log ${input.prefix}_dupmetrics.log --noPhoneHome
        ""","MarkDups2"
    }
}
