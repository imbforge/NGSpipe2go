VariantScoreRecalibration = {
    doc title: "GATK Variant Quality Score Recalibration",
        desc:  "Calculate VQSLOD scores for further filtering variants. See https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = VariantScoreRecalibration_vars.outdir

    def VariantRecalibrator_INDEL_FLAGS =
        " -mode INDEL --trustAllPolymorphic -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP " +
        (VariantScoreRecalibration_vars.bwa_ref        ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.mills_variants ? " -resource:mills,known=false,training=true,truth=true,prior=12 " + VariantScoreRecalibration_vars.mills_variants : "" ) +
        (VariantScoreRecalibration_vars.known_variants ? " -resource:dbsnp,known=true,training=false,truth=false,prior=2 " + VariantScoreRecalibration_vars.known_variants : "" ) +
        (VariantScoreRecalibration_vars.max_gaussians_indels ? " --maxGaussians " + VariantScoreRecalibration_vars.max_gaussians_indels : "" )

    def ApplyRecalibration_INDEL_FLAGS =
        " -mode INDEL " +
        (VariantScoreRecalibration_vars.bwa_ref        ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.known_variants ? " -ts_filter_level " + VariantScoreRecalibration_vars.indel_filter_level : "" )

    def VariantRecalibrator_SNP_FLAGS =
        " -mode SNP --trustAllPolymorphic -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP " +
        (VariantScoreRecalibration_vars.bwa_ref                   ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.hapmap_variants           ? " -resource:hapmap,known=false,training=true,truth=true,prior=15 " + VariantScoreRecalibration_vars.hapmap_variants           : "" ) +
        (VariantScoreRecalibration_vars.omni_variants             ? " -resource:omni,known=false,training=true,truth=true,prior=12 "   + VariantScoreRecalibration_vars.omni_variants             : "" ) +
        (VariantScoreRecalibration_vars.thousand_genomes_variants ? " -resource:1000G,known=false,training=true,truth=false,prior=10 " + VariantScoreRecalibration_vars.thousand_genomes_variants : "" ) +
        (VariantScoreRecalibration_vars.known_variants            ? " -resource:dbsnp,known=true,training=false,truth=false,prior=7 "  + VariantScoreRecalibration_vars.known_variants            : "" ) +
        (VariantScoreRecalibration_vars.max_gaussians_indels      ? " --maxGaussians " + VariantScoreRecalibration_vars.max_gaussians_snps : "" )

    def ApplyRecalibration_SNP_FLAGS =
        " -mode SNP " +
        (VariantScoreRecalibration_vars.bwa_ref        ? " -R " + VariantScoreRecalibration_vars.bwa_ref : "" ) +
        (VariantScoreRecalibration_vars.known_variants ? " --ts_filter_level " + VariantScoreRecalibration_vars.indel_filter_level : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform (".vcf.gz") to (".indels.recal", ".indels.tranches", ".snps.recal", ".snps.tranches", ".vqsr.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${VariantScoreRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T VariantRecalibrator --input $input --recal_file $output1 --tranches_file $output2 $VariantRecalibrator_INDEL_FLAGS &&
            java ${VariantScoreRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T ApplyRecalibration --input $input --recal_file $output1 --tranches_file $output2 -o \${TMP}/\$(basename ${input}).indel.vqsr.vcf.gz $ApplyRecalibration_INDEL_FLAGS &&
            java ${VariantScoreRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T VariantRecalibrator --input \${TMP}/\$(basename ${input}).indel.vqsr.vcf.gz --recal_file $output3 --tranches_file $output4 $VariantRecalibrator_SNP_FLAGS &&
            java ${VariantScoreRecalibration_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T ApplyRecalibration --input \${TMP}/\$(basename ${input}).indel.vqsr.vcf.gz --recal_file $output3 --tranches_file $output4 -o $output5 $ApplyRecalibration_SNP_FLAGS
        ""","VariantScoreRecalibration"
    }
    
}
