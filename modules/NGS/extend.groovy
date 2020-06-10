extend = {
    doc title: "extend",
        desc:  "Extend read length to the average fragment size",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir=extend_vars.outdir

    def SAMTOOLS_SORT_FLAGS = "-O bam -@ " + extend_vars.samtools_threads

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("extend")

    transform(".bam") to ("_ext.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            CHRSIZES="\${TMP}/\$(basename ${input.prefix}).extend.chrsizes"  &&
            samtools idxstats ${input} | cut -f1-2 > "\${CHRSIZES}" &&
            bedtools bamtobed -split -i $input | bedtools slop -g "\${CHRSIZES}" -l 0 -r ${extend_vars.fraglen} -s | bedtools bedtobam -ubam -g "\${CHRSIZES}" | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix) - > $output &&
            samtools index $output
        ""","extend"
    }
}

