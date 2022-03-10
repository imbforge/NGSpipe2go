VariantScoreRecalibration = {
    doc title: "GATK Variant Quality Score Recalibration",
        desc:  "Filtering technique applied on the variant callset that uses machine learning to identify annotation profiles of variants that are likely to be real and assigns a VQSLOD score to each variant that is much more reliable than the QUAL score calculated by the caller. The obtained model based on training variants is then applied to the data to filter out probable artifacts from the callset. See https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-",
        constraints: "Comprehensive reference data panels required. There are many things that can cause this tool to fail and it rarely provides a useful error message.",
        bpipe_version: "tested with bpipe 0.9.9.8.slurm",
        author: "Sergi Sayols, modified by Frank Rühle"

    output.dir = VariantScoreRecalibration_vars.outdir

    def VariantRecalibrator_INDEL_FLAGS =
        " -mode INDEL --trust-all-polymorphic -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP " +
        (VariantScoreRecalibration_vars.bwa_ref        ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.mills_variants ? " -resource:mills,known=false,training=true,truth=true,prior=12 " + VariantScoreRecalibration_vars.mills_variants : "" ) +
        (VariantScoreRecalibration_vars.known_variants ? " -resource:dbsnp,known=true,training=false,truth=false,prior=2 " + VariantScoreRecalibration_vars.known_variants : "" ) +
        (VariantScoreRecalibration_vars.max_gaussians_indels ? " --max-gaussians " + VariantScoreRecalibration_vars.max_gaussians_indels : "" )

    def ApplyVQSR_INDEL_FLAGS =
        " -mode INDEL " +
        (VariantScoreRecalibration_vars.bwa_ref        ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.known_variants ? " --truth-sensitivity-filter-level " + VariantScoreRecalibration_vars.indel_filter_level : "" )

    def VariantRecalibrator_SNP_FLAGS =
        " -mode SNP --trust-all-polymorphic -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP " +
        (VariantScoreRecalibration_vars.bwa_ref                   ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.hapmap_variants           ? " -resource:hapmap,known=false,training=true,truth=true,prior=15 " + VariantScoreRecalibration_vars.hapmap_variants           : "" ) +
        (VariantScoreRecalibration_vars.omni_variants             ? " -resource:omni,known=false,training=true,truth=true,prior=12 "   + VariantScoreRecalibration_vars.omni_variants             : "" ) +
        (VariantScoreRecalibration_vars.thousand_genomes_variants ? " -resource:1000G,known=false,training=true,truth=false,prior=10 " + VariantScoreRecalibration_vars.thousand_genomes_variants : "" ) +
        (VariantScoreRecalibration_vars.known_variants            ? " -resource:dbsnp,known=true,training=false,truth=false,prior=7 "  + VariantScoreRecalibration_vars.known_variants            : "" ) +
        (VariantScoreRecalibration_vars.max_gaussians_indels      ? " --max-gaussians " + VariantScoreRecalibration_vars.max_gaussians_snps : "" )

    def ApplyVQSR_SNP_FLAGS =
        " -mode SNP " +
        (VariantScoreRecalibration_vars.bwa_ref        ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.known_variants ? " --truth-sensitivity-filter-level " + VariantScoreRecalibration_vars.indel_filter_level : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    transform (".vcf.gz") to (".indels.recal", ".indels.tranches", ".indels.plots.R", ".snps.recal", ".snps.tranches", ".snps.plots.R", ".vqsr.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${VariantScoreRecalibration_vars.java_flags}" VariantRecalibrator -V $input -O $output1 --tranches-file $output2 --rscript-file $output3 $VariantRecalibrator_INDEL_FLAGS &&
            gatk --java-options "${VariantScoreRecalibration_vars.java_flags}" ApplyVQSR -V $input --recal-file $output1 --tranches-file $output2 -O \${TMP}/\$(basename ${input}).indel.vqsr.vcf.gz $ApplyVQSR_INDEL_FLAGS &&
            gatk --java-options "${VariantScoreRecalibration_vars.java_flags}" VariantRecalibrator -V \${TMP}/\$(basename ${input}).indel.vqsr.vcf.gz -O $output4 --tranches-file $output5 --rscript-file $output6 $VariantRecalibrator_SNP_FLAGS &&
            gatk --java-options "${VariantScoreRecalibration_vars.java_flags}" ApplyVQSR -V \${TMP}/\$(basename ${input}).indel.vqsr.vcf.gz --recal-file $output4 --tranches-file $output5 -O $output7 $ApplyVQSR_SNP_FLAGS
        ""","VariantScoreRecalibration"
    }
    
}
