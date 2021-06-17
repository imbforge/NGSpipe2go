Bowtie_se = {
    doc title: "Bowtie SE alignment",
        desc:  "Align single end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        author: "Sergi Sayols, Anke Busch"

    output.dir = Bowtie_se_vars.mapped
    def BOWTIE_FLAGS =
        " -q --sam" +
        (Bowtie_se_vars.quals       ? " "    + Bowtie_se_vars.quals       : "") +
        (Bowtie_se_vars.best        ? " --best --strata --tryhard --chunkmbs 256" : "") +
        (Bowtie_se_vars.threads     ? " -p " + Bowtie_se_vars.threads     : "") +
        (Bowtie_se_vars.mm          ? " -v " + Bowtie_se_vars.mm          : "") +
        (Bowtie_se_vars.multireport ? " -M " + Bowtie_se_vars.multireport : "") +
        (Bowtie_se_vars.extra   ? " "    + Bowtie_se_vars.extra    : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (Bowtie_se_vars.samtools_threads ? " -@ " + Bowtie_se_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(module:"Bowtie_se", branch:branch, branch_outdir:"")

    transform(".fastq.gz") to (".bam") {
        def SAMPLENAME = output.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&

            echo 'BOWTIE_FLAGS' $BOWTIE_FLAGS > $output.dir/\${SAMPLENAME_BASE}.bowtie.log &&
            echo 'BOWTIE_REF' $Bowtie_se_vars.ref >> $output.dir/\${SAMPLENAME_BASE}.bowtie.log && 

            zcat $input | bowtie $BOWTIE_FLAGS $Bowtie_se_vars.ref - 2>> $output.dir/\${SAMPLENAME_BASE}.bowtie.log | awk '{if (\$1~/^@/) print; else {if(\$5 == 255) print \$0"\tNH:i:1"; else print \$0"\tNH:i:2";}}' | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix)_bowtie1_sort - -o $output
           ""","Bowtie_se"
    }
}
