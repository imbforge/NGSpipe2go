// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/merge_SJ_tab.vars.groovy"

FilterAndMergeSJtab = {
    doc title: "FilterAndMergeSJtab",
        desc: "GATK variant calling suggests 2-step STAR mapping for RNA-seq. In this steps all splice junctions files are collected, filtered and merged. Based on https://code.google.com/p/rna-star/issues/detail?id=7",
        constraints: "STAR STAR_2.4.2a",
        author: "Antonio Domingues"

    output.dir = OUTDIR_2ND_INDEX

    def PREAMBLE = get_preamble("FilterAndMergeSJtab")

    produce("SJ.out.tab.Pass1.sjdb"){
        exec """
            ${PREAMBLE} &&

            cat $inputs | awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' > ${OUTDIR_2ND_INDEX}/SJ.out.tab.Pass1.sjdb
        ""","FilterAndMergeSJtab"
    }
}
