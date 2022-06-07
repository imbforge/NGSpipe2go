count_breaks = {
    doc title: "Count breaks per position",
    desc:  "Count breals. Executes a special UMI filtering step for breaktag (though UMIs are not mandatory)",
    constraints: "none",
    author: "Sergi Sayols"

    output.dir = count_breaks_vars.outdir

    def samtools_view_FLAGS = "-b" +
        (count_breaks_vars.paired ? " -f 0x0040" : "")

    def count_breaks_FLAGS =
        (count_breaks_vars.threads ? " -t " + count_breaks_vars.threads : "" )

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"]) + " && " +
                   prepare_tool_env("python"  , tools["python"]["version"]  , tools["python"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bam") to(".bed.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            tmpfile=\$(mktemp -p ${TMP}) &&

            samtools view $samtools_view_FLAGS $input | \
              bedtools bamtobed -i - | \
              awk 'BEGIN{OFS=","} {sub(/^.+_/,"",\$4); if(\$6 == "+") { print \$1,\$2-1,\$2,\$6,\$4 } else { print \$1,\$3,\$3+1,\$6,\$4 }}' | \
              sort -t, --parallel=$count_breaks_vars.threads -k1,1 -k2,2g -k3,3g | \
              uniq -c | \
              awk 'BEGIN{OFS=","} { print \$2,\$1 }' > \$tmpfile &&

            python3 ${PIPELINE_ROOT}/tools/breaktag/umi_filtering.py \$tmpfile | \
              cut -f-4 | \
              awk '{ if (\$4 == "-") {\$2=\$2-1;\$3=\$3-1}; print }' | \
              sort --parallel=$count_breaks_vars.threads -k1,1 -k2,2g -k3,3g | \
              uniq -c | \
              awk 'BEGIN{OFS="\t"} { print \$2,\$3,\$4,".",\$1,\$5 }' | \
              gzip -c > $output
            ""","count_breaks"
    }
}


