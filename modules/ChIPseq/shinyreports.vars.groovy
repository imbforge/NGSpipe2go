shinyReports_vars=[
    project         : PROJECT,              //project directory
    log             : LOGS,                 //where the logs lie
    qc              : QC,                   //where the QC lie
    res             : RESULTS,              //where the results lie
    prefix          : ESSENTIAL_SAMPLE_PREFIX,
    targets         : "targets.txt",
    bowtie_log      : LOGS + "/bowtie_se",  //where the Bowtie logs lie
    bamindex_log    : LOGS + "/BAMindexer", //where the Samtools/BamIndexer logs lie
    markdups_log    : LOGS + "/MarkDups",   //where the MarkDups logs lie
    extend_log      : LOGS + "/extend",     //where the extend/BedTools logs lie
    fastqc          : FastQC_vars.outdir,   //where the Fastqc logs lie
    fastqc_log      : LOGS + "/FastQC",     //where the Fastqc logs lie
    ipstrength      : ipstrength_vars.outdir,   //where the IPstrength files lie
    ipstrength_log  : LOGS + "/ipstrength",     //where the IPstrength/R logs lie
    pbc             : pbc_vars.outdir,          //where the PBC files lie
    phantompeak     : phantompeak_vars.outdir,  //where the PhantomPeak files lie
    phantom_log     : LOGS + "/phantompeak",    //where the PhantomPeak/R logs lie
    bustard         : QC + "/DemultiplexedBustardSummary.xml",    //where the bustard xml file lies
    macs2           : macs2_vars.outdir,    //where the MACS2 results lie
    macs2_log       : LOGS + "/macs2",      //where the macs2 logs lie
    blacklist_filter: blacklist_filter_vars.outdir + "/peaks_detected_table.csv",
    plots_column    : "4",                  //number of columns to splits the plot grids (ipstrength, phantompeaks...). Min=2
    peak_annotation : peak_annotation_vars.outdir, // where the peak annotation results lie
    great           : GREAT_vars.outdir,    // where the GO enrichment results lie
    trackhub_done   : PROJECT + "/trackhub.done",   // contains trackhub URL
    tool_versions   : PIPELINE_ROOT + "/pipelines/ChIPseq/tools.groovy" //where the tool versions listed
]
