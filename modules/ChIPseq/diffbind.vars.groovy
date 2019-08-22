diffbind_vars=[
    targets          : "targets_diffbind.txt",   //targets file. Check the bin directory for the format
    contrasts        : "contrasts_diffbind.txt", //contrasts file. Check the bin directory for the format
    cwd              : PROJECT,                  //current working directory
    outdir           : RESULTS + "/diffbind",    //output directory
    bams             : MAPPED,                   // directory with the bam files
    fragsize         : Integer.toString(ESSENTIAL_FRAGLEN), // average fragment size
    substractcontrol : true,                     // substract input
    fulllibrarysize  : true,                     // use total number of reads in bam for normalization (FALSE=only peaks)
    tagwisedispersion: true,                     // calculate dispersion tagwise (use FALSE if no replicates)
    annotate         : true,                     // annotate peaks after differential bindign analysis?
    tss              : "'c(-3000,3000)'",        // region around the TSS to be considered as promoter
    txdb             : ESSENTIAL_TXDB,           // Bioconductor transcript database, for annotation
    annodb           : ESSENTIAL_ANNODB,         // Bioconductor gene annotation database
    paired           : (ESSENTIAL_PAIRED == "yes"), // Paired end experiment?
    extra            : ""                        //extra parms to sent to the tool
]
