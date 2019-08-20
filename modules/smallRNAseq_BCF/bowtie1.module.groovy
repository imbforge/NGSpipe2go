// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/bowtie1.vars.groovy"

Bowtie_se = {
    doc title: "Bowtie SE alignment",
        desc:  "Align single end reads",
        constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
        author: "Sergi Sayols, Anke Busch"

    output.dir = MAPPED
    def BOWTIE_FLAGS = " -q --sam" +
                       " " + BOWTIE_QUALS +
                       " " + BOWTIE_BEST +
                       " -p " + Integer.toString(BOWTIE_THREADS) +
                       " -v " + Integer.toString(BOWTIE_MM) +
                       " -M " + Integer.toString(BOWTIE_MULTIREPORT) 

    def TOOL_ENV = prepare_tool_env("bowtie", tools["bowtie"]["version"], tools["bowtie"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("Bowtie_se")

    transform(".fastq.gz") to (".bam") {
        def SAMPLENAME = output.prefix
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&

            echo 'BOWTIE_FLAGS' $BOWTIE_FLAGS > $output.dir/\${SAMPLENAME_BASE}.bowtie.log &&
            echo 'BOWTIE_REF' $BOWTIE_REF >> $output.dir/\${SAMPLENAME_BASE}.bowtie.log && 

            zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - 2>> $output.dir/\${SAMPLENAME_BASE}.bowtie.log | awk '{if (\$1~/^@/) print; else {if(\$5 == 255) print \$0"\tNH:i:1"; else print \$0"\tNH:i:2";}}' | samtools view -bhSu - | samtools sort -@ $BOWTIE_THREADS -T \${TMP}/\$(basename $output.prefix)_bowtie1_sort - -o $output
           ""","Bowtie_se"
    }
}
