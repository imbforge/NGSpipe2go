shinyReports_vars=[
    project         : PROJECT,          //project directory
    expdesign       : ESSENTIAL_EXPDESIGN,
    log             : LOGS,             //where the logs lie
    qc              : QC,               //where the QC lie
    res             : RESULTS,          //where the results lie
    fastqc_out      : FastQC_vars.outdir,      //where the Fastqc output lie
    fastqc_log      : LOGS + "/FastQC",        //where the Fastqc logs lie
    bamindex_log    : LOGS + "/BAMindexer",    //where the Samtools/BamIndexer logs lie
    bustard         : QC + "/DemultiplexedBustardSummary.xml",  //where the bustard xml file lies
    bam2bw_log      : LOGS + "/bam2bw",              //where the Bam2BW logs lie
    pear_out        :  pear_vars.outdir,
    pear_log        :  pear_vars.logdir,
    umiextract_out  :  AddUMIBarcodeToFastq_vars.outdir,
    umiextract_log  :  AddUMIBarcodeToFastq_vars.logdir,
    barcode_count_out   :  barcode_count_vars.outdir,
    barcode_count_log   :  barcode_count_vars.logdir,
    plots_column    : 3,                             //number of columns to splits the plot grids (dupradar, genebodycov...). Min=2L. L=integer in R
    trackhub_done   : PROJECT + "/trackhub.done",    //contains trackhub URL
    tool_versions   : collectToolVersions_vars.outdir + "/tool_versions.txt", //where the tool versions listed
    target          : new File(PIPELINE_ROOT + "/pipelines/DNAampliconseq/targets.txt").getCanonicalPath()
]
