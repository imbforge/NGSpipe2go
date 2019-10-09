// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load  PIPELINE_ROOT + "/modules/DNAseq/variantcallHC.vars.groovy"

VariantCallHC = {
    doc title: "GATK Variant Calling HC",
        desc:  "Call variants in BAM files using GATK HaplotypeCaller",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Oliver Drechsel"

    output.dir = VariantCallHC_vars.outdir

    def HaplotypeCaller_FLAGS =
        (VariantCallHC_vars.call_region    ? " -L "      + VariantCallHC_vars.call_region    : "" ) +
        (VariantCallHC_vars.threads        ? " -nct "    + VariantCallHC_vars.threads        : "" ) +
        (VariantCallHC_vars.bwa_ref        ? " -R "      + VariantCallHC_vars.bwa_ref        : "" ) +
        (VariantCallHC_vars.known_variants ? " --dbsnp " + VariantCallHC_vars.known_variants : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("VariantCallHC")

    transform (".dupmarked.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${VariantCallHC_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T HaplotypeCaller $HaplotypeCaller_FLAGS -I $input -o $output
        ""","VariantCallHC"
    }
}

