VariantFiltration = {
    doc title: "GATK VariantFiltration",
        desc:  "Hard filter a cohort callset with GATK VariantFiltration. Filter criteria are applied separately to SNPs and Indels (mixed-type variants are treated like Indels) and resulting vcf files are sorted and merged afterwards. Hard-filtering is useful when the data cannot support VQSR or when an analysis requires manual filtering (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering).",
        constraints: "",
        bpipe_version: "tested with 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = VariantFiltration_vars.outdir

    def VariantFiltration_INDEL_FLAGS =
        (VariantFiltration_vars.bwa_ref           ? " -R " + VariantFiltration_vars.bwa_ref : "" ) +
        (VariantFiltration_vars.indel_filter_QD   ? " -filter \"" + VariantFiltration_vars.indel_filter_QD + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_QD.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.indel_filter_QUAL ? " -filter \"" + VariantFiltration_vars.indel_filter_QUAL + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_QUAL.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.indel_filter_SOR  ? " -filter \"" + VariantFiltration_vars.indel_filter_SOR + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_SOR.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.indel_filter_FS   ? " -filter \"" + VariantFiltration_vars.indel_filter_FS + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_FS.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.indel_filter_MQ   ? " -filter \"" + VariantFiltration_vars.indel_filter_MQ + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_MQ.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.indel_filter_MQRankSum      ? " -filter \"" + VariantFiltration_vars.indel_filter_MQRankSum + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_MQRankSum.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.indel_filter_ReadPosRankSum ? " -filter \"" + VariantFiltration_vars.indel_filter_ReadPosRankSum + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.indel_filter_ReadPosRankSum.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) 

    def VariantFiltration_SNP_FLAGS =
        (VariantFiltration_vars.bwa_ref         ? " -R " + VariantFiltration_vars.bwa_ref : "" ) +
        (VariantFiltration_vars.snp_filter_QD   ? " -filter \"" + VariantFiltration_vars.snp_filter_QD + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_QD.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.snp_filter_QUAL ? " -filter \"" + VariantFiltration_vars.snp_filter_QUAL + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_QUAL.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.snp_filter_SOR  ? " -filter \"" + VariantFiltration_vars.snp_filter_SOR + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_SOR.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.snp_filter_FS   ? " -filter \"" + VariantFiltration_vars.snp_filter_FS + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_FS.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.snp_filter_MQ   ? " -filter \"" + VariantFiltration_vars.snp_filter_MQ + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_MQ.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.snp_filter_MQRankSum      ? " -filter \"" + VariantFiltration_vars.snp_filter_MQRankSum + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_MQRankSum.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) +
        (VariantFiltration_vars.snp_filter_ReadPosRankSum ? " -filter \"" + VariantFiltration_vars.snp_filter_ReadPosRankSum + "\"" + 
        " --filter-name \"" + VariantFiltration_vars.snp_filter_ReadPosRankSum.replaceAll("[^a-zA-Z0-9.-]","") + "\"" : "" ) 

    def VariantsToTable_FLAGS = "--show-filtered -F CHROM -F POS -F ID -F REF -F ALT -F FILTER -F TYPE -F EVENTLENGTH -F QUAL -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F MULTI-ALLELIC -F NSAMPLES -F NCALLED -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum"

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"]) + " && " +
                   prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    transform (".vcf.gz") to (".filtered.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            gatk --java-options "${VariantFiltration_vars.java_flags}" SelectVariants -V $input -select-type INDEL -select-type MIXED -O \${TMP}/indel_\$(basename ${input}) &&
            gatk --java-options "${VariantFiltration_vars.java_flags}" VariantFiltration -V \${TMP}/indel_\$(basename ${input}) -O \${TMP}/indel_filtered_\$(basename ${input}) $VariantFiltration_INDEL_FLAGS &&
            gatk --java-options "${VariantFiltration_vars.java_flags}" SelectVariants -V $input -select-type SNP -O \${TMP}/snp_\$(basename ${input}) &&
            gatk --java-options "${VariantFiltration_vars.java_flags}" VariantFiltration -V \${TMP}/snp_\$(basename ${input}) -O \${TMP}/snp_filtered_\$(basename ${input}) $VariantFiltration_SNP_FLAGS &&
            java ${VariantFiltration_vars.java_flags} -jar \${PICARD} SortVcf I=\${TMP}/indel_filtered_\$(basename ${input}) O=\${TMP}/indel_filtered_sorted_\$(basename ${input}) &&
            java ${VariantFiltration_vars.java_flags} -jar \${PICARD} SortVcf I=\${TMP}/snp_filtered_\$(basename ${input}) O=\${TMP}/snp_filtered_sorted_\$(basename ${input}) &&
            java ${VariantFiltration_vars.java_flags} -jar \${PICARD} MergeVcfs I=\${TMP}/indel_filtered_sorted_\$(basename ${input}) I=\${TMP}/snp_filtered_sorted_\$(basename ${input}) O=$output &&
            gatk --java-options "${VariantFiltration_vars.java_flags}" VariantsToTable -V $output -O ${output.dir}/filteredVariants.table $VariantsToTable_FLAGS
            
        ""","VariantFiltration"
    }
    
}
