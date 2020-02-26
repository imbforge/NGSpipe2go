// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load  PIPELINE_ROOT + "/modules/DNAseq/realignment.vars.groovy"

IndelRealignment = {
    doc title: "GATK IndelRealignment",
        desc:  "Realign BAM files at Indel positions, using GATK",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = IndelRealignment_vars.outdir

    def RealignerTargetCreator_FLAGS =
        (IndelRealignment_vars.threads        ? " -nt "     + IndelRealignment_vars.threads : "" ) +
        (IndelRealignment_vars.bwa_ref        ? " -R "      + IndelRealignment_vars.bwa_ref : "" ) +
        (IndelRealignment_vars.mills_variants ? " -known " + IndelRealignment_vars.mills_variants : "" )

    def IndelRealignment_FLAGS =
        (IndelRealignment_vars.bwa_ref        ? " -R "      + IndelRealignment_vars.bwa_ref : "" ) +
        (IndelRealignment_vars.mills_variants ? " -known " + IndelRealignment_vars.mills_variants : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("IndelRealignment")

    transform (".bam") to (".realignment.targets.bed", ".realigned.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${IndelRealignment_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T RealignerTargetCreator $RealignerTargetCreator_FLAGS -I $input -o $output1 &&
            java ${IndelRealignment_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T IndelRealigner $IndelRealignment_FLAGS -I $input -targetIntervals $output1 -o $output2
       ""","IndelRealignment"
    }
    
}
