count_breaks_strandless = {
    doc title: "Count breaks per position",
    desc:  "Count breaks regardless of the strand where the read points them to be",
    constraints: "Expect to have perl installed",
    author: "Sergi Sayols"

    output.dir = count_breaks_strandless_vars.outdir

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bed.gz") to(".stradnless.bed.gz") {
        exec """
            ${PREAMBLE} &&

            zcat $input | \
              perl -aln -e 'if(\$F[0]==\$F_prev[0] && \$F[1]==\$F_prev[1] && \$F[2]==\$F_prev[2]){ \$F_prev[5]+=\$F[5]; } else { \$F[5]="*"; print join("\t", @F_prev); @F_prev=@F; } END{ print join("\t", @F) }' | \
              tail -n +2 | \
              gzip -c > $output
        ""","count_breaks_strandless"
    }
}


