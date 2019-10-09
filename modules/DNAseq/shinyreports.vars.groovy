shinyReports_vars=[
    project               : PROJECT,                 //project directory
    log                   : LOGS,                    //where the logs lie
    qc                    : QC,                      //where the QC lie
    res                   : RESULTS,                 //where the results lie
    prefix                : ESSENTIAL_SAMPLE_PREFIX, //standard sample prefix
    fastqc_out            : FastQC_vars.outdir,      //where the Fastqc output lie
    fastqc_log            : LOGS + "/FastQC",        //where the Fastqc logs lie
    bwa_log               : LOGS + "/BWA_pe",        //where the BWA (samtools flagstat) logs lie
    bwa_suffix            : ".bam.log",              //extension given to the BWA log files
    markdups_log          : LOGS + "/MarkDups",      //where the picard MarkDuplicates logs lie
    gatkug_log            : LOGS + "/VariantCallUG", //where the GATK UnifiedGenotyper logs lie
    gatkug_suffix         : ".UG.vcf.gz.log",        //extension given to the GATK UnifiedGenotyper log files
    gatkhc_log            : LOGS + "/VariantCallHC", //where the GATK HaplotypeCaller logs lie
    gatkhc_suffix         : ".HC.vcf.gz.log",        //extension given to the GATK Haplotypecaller log files
    gatkvarianteval       : VariantEval_vars.outdir, //location of GATK variantEval results
    gatkvarianteval_suffix: ".report",               //extension given to the GATK VariantEval output files
    res_gatkhc_suffix     : ".HC.vcf.gz",            //extension of the GATK HaplotypeCaller output files
    tool_versions         : PIPELINE_ROOT + "/pipelines/ChIPseq/tools.groovy" //where the tool versions listed
]
