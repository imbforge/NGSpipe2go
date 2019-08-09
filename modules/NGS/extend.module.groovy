// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/extend.vars.groovy"

extend = {
    doc title: "extend",
        desc:  "Extend read length to the average fragment size",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir=MAPPED

    def SAMTOOLS_SORT_FLAGS = "-O bam " + EXTEND_SAMTOOLS_THREADS

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])

    transform(".bam") to ("_ext.bam") {
        exec """
            ${TOOL_ENV} &&

            if [ ! -d $TMP ]; then
                mkdir -p $TMP;
            fi &&

            CHRSIZES=${TMP}/\$(basename ${input.prefix}).extend.chrsizes  &&
            samtools idxstats ${input} | cut -f1-2 > \${CHRSIZES} &&
            bedtools bamtobed -split -i $input |
            bedtools slop -g \$CHRSIZES -l 0 -r $EXTEND_FRAGLEN -s |
            bedtools bedtobam -ubam -g \$CHRSIZES |
            samtools sort $SAMTOOLS_SORT_FLAGS -T $TMP/\$(basename $output.prefix) - > $output &&
            samtools index $output
        ""","extend"
    }
}

