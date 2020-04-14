shinyReports_vars=[
    project         : PROJECT,          //project directory
    expdesign       : ESSENTIAL_EXPDESIGN,
    target          : new File(PIPELINE_ROOT + "/pipelines/DNAampliconseq/targets.txt").getCanonicalPath(),
    log             : LOGS,             //where the logs lie
    qc              : QC,               //where the QC lie
    res             : RESULTS,          //where the results lie
    fastqc_out      : FastQC_vars.outdir,      //where the Fastqc output lie
    fastqc_log      : LOGS + "/FastQC",        //where the Fastqc logs lie
    patternFastQC   : "", // pattern for sample names to select do display in report
    patternFastQC2  : "cutadapt", // pattern2 for sample names to select do display in report (e.g. after adapter trimming)
    patternUmitools : "nested", // pattern2 for sample names to select do display in report (e.g. after adapter trimming)
    maxno           : "Inf", //the maximum plot number will be restricted accordingly to the first 'maxno' plots.
    colorByFactor   : "sub_experiment", // default variables for grouping and plotting (max 2)
    bamindex_log    : LOGS + "/BAMindexer",    //where the Samtools/BamIndexer logs lie
    bustard         : QC + "/DemultiplexedBustardSummary.xml",  //where the bustard xml file lies
    bam2bw_log      : LOGS + "/bam2bw",              //where the Bam2BW logs lie
    pear_out        :  pear_vars.outdir,
    pear_log        :  pear_vars.logdir,
    umiextract_out  :  AddUMIBarcodeToFastq_vars.outdir,
    umiextract_log  :  AddUMIBarcodeToFastq_vars.logdir,
    umiextract_logWL    :  AddUMIBarcodeToFastq_vars.logdirWL,
    barcode_count_out   :  barcode_count_vars.outdir,
    barcode_count_log   :  barcode_count_vars.logdir,
    MPSprofiling_out    :  MPSprofiling_vars.outdir,
    MPSprofiling_log    :  MPSprofiling_vars.logdir,
    MPS_threshold_rel_countssum      : MPSprofiling_vars.threshold_rel_countssum,
    MPS_excludeSeqsNotInAllFractions : MPSprofiling_vars.excludeSeqsNotInAllFractions,
    maxTableRows    : "10000",  // maximum row number to be displayed in PSI result tables. "Inf" mean no restriction.
    plots_column    : "2",                             //number of columns to splits the plot grids (dupradar, genebodycov...). Min=2L. L=integer in R
    trackhub_done   : PROJECT + "/trackhub.done",    //contains trackhub URL
    tool_versions   : collectToolVersions_vars.outdir + "/tool_versions.txt" //where the tool versions listed
]
