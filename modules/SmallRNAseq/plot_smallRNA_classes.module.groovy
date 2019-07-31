load MODULE_FOLDER + "SmallRNAseq/plot_smallRNA_classes.vars.groovy"

PlotSmallRNAclasses = {
    doc title: "PlotSmallRNAclasses",
    desc: "Plots the number of reads in each filtered file, normalized to mapped reads.",
    constraints: "Reuires a table with library sizes.",
    author: "António Domingues"

    output.dir = PLOT_CLASSES_OUTDIR

    produce("smallRNAClassesCount.txt"){
        exec """
            module load R/${R_VERSION} &&

            bams=`find ${PLOT_CLASSES_OUTDIR} -name "*.bam"` &&

            Rscript $PLOT_SMALL_RNA_TOOL_PATH --bams $bams --libsizes ${PLOT_CLASSES_MAPPED} --out $output.dir 

      ""","PlotSmallRNAclasses"
    }
}
