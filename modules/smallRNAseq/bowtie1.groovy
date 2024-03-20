bowtie1_sRNA = {
    doc title: "Bowtie1 smallRNA alignment",
        desc:  "Align single end smallRNA reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        author: "Sergi Sayols, Anke Busch"

    output.dir = bowtie1_sRNA_vars.mapped

    // create the log folder if it doesn't exist
    def Bowtie_LOGDIR = new File(bowtie1_sRNA_vars.logdir)
    if (!Bowtie_LOGDIR.exists()) {
        Bowtie_LOGDIR.mkdirs()
    }

    def BOWTIE_FLAGS =
        " -q --sam" +
        (bowtie1_sRNA_vars.quals       ? " "    + bowtie1_sRNA_vars.quals       : "") +
        (bowtie1_sRNA_vars.best        ? " --best --strata --tryhard --chunkmbs 256" : "") +
        (bowtie1_sRNA_vars.threads     ? " -p " + bowtie1_sRNA_vars.threads     : "") +
        (bowtie1_sRNA_vars.mm          ? " -v " + bowtie1_sRNA_vars.mm          : "") +
        (bowtie1_sRNA_vars.multireport ? " -M " + bowtie1_sRNA_vars.multireport : "") +
        (bowtie1_sRNA_vars.extra       ? " "    + bowtie1_sRNA_vars.extra       : "")

    def SAMTOOLS_VIEW_FLAGS = "-bhSu "
    def SAMTOOLS_SORT_FLAGS = "-O bam " +
        (bowtie1_sRNA_vars.samtools_threads ? " -@ " + bowtie1_sRNA_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    transform(".fastq.gz") to (".bam") {
        def SAMPLENAME = output.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&

            echo 'BOWTIE_FLAGS' $BOWTIE_FLAGS > $Bowtie_LOGDIR/\${SAMPLENAME_BASE}.bowtie.log &&
            echo 'BOWTIE_REF' $bowtie1_sRNA_vars.ref >> $Bowtie_LOGDIR/\${SAMPLENAME_BASE}.bowtie.log && 

            zcat $input | bowtie $BOWTIE_FLAGS $bowtie1_sRNA_vars.ref - 2>> $Bowtie_LOGDIR/\${SAMPLENAME_BASE}.bowtie.log | awk '{if (\$1~/^@/) print; else {if(\$5 == 255) print \$0"\tNH:i:1"; else print \$0"\tNH:i:2";}}' | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/\$(basename $output.prefix)_bowtie1_sort - -o $output
           ""","bowtie1_sRNA"
    }
}
