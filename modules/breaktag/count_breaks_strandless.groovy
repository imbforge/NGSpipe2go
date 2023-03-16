count_breaks_strandless = {
    doc title: "Count breaks per position",
    desc:  "Count breaks regardless of the strand where the read points them to be",
    constraints: "Expect to have perl installed",
    author: "Sergi Sayols"

    output.dir = count_breaks_strandless_vars.outdir

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bed.gz") to(".strandless.bed.gz") {
        exec """
            ${PREAMBLE} &&

            zcat $input | \
              perl -aln -e 'if(\$F[0]==\$F0[0] && \$F[1]==\$F0[1] && \$F[2]==\$F0[2]){ \$F0[4]+=\$F[4]; } else { \$F0[5]="*"; print join("\t", @F0); @F0=@F; } END{ \$F[5]="*"; print join("\t", @F) }' | \
              tail -n +2 | \
              gzip -c > $output
        ""","count_breaks_strandless"
    }
}


