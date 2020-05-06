DE_DESeq2_MM_vars=[
    targets  : DE_DESeq2_vars.targets,        //targets file. Check the bin directory for the format
    contrasts: DE_DESeq2_vars.contrasts,      //contrasts file. Check the bin directory for the format
    mmatrix  : DE_DESeq2_vars.mmatrix,        //formula for the linear model
    filter   : DE_DESeq2_vars.filter,         //whether to perform automatic independent filtering of lowly expressed genes
    prefix   : DE_DESeq2_vars.prefix,         //prefix to be removed from file names
    suffix   : DE_DESeq2_vars.suffix,         //suffix to be removed from file names
    dupradar_outdir: dupRadar_vars.outdir,    //where the dupradar left the output tables
    cwd      : RESULTS + "/subread-count_MM", //directory where the .tsv files with the gene counts are located
    outdir   : RESULTS + "/DE_DESeq2_MM",     //output filename base pattern. If you change it here, change it also in the module file
    genes    : DE_DESeq2_vars.genes,          //gtf file to calculate RPKM
    pattern  : "\"\\.readcounts.tsv\"",       //pattern for filtering count files in cwd
    extra    : DE_DESeq2_vars.extra           //extra parms to sent to the tool
]
