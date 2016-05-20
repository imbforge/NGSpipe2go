//  rule to filter and merge SJ.out.tabs to prepare 2-step star mapping, version 1
FilterAndMergeSJtab = {
   doc title: "FilterAndMergeSJtab",
   desc: "GATK variant calling suggests 2-step STAR mapping for RNA-seq. In this steps all splice junctions files are collected, filtered and merged. Based on https://code.google.com/p/rna-star/issues/detail?id=7",
   constraints: "STAR STAR_2.4.2a",
   author: "Antonio Domingues"

   output.dir = OUTDIR_2ND_INDEX

   produce("SJ.out.tab.Pass1.sjdb"){
      exec """
                        cat $inputs | awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' > ${OUTDIR_2ND_INDEX}/SJ.out.tab.Pass1.sjdb
      ""","FilterAndMergeSJtab"
   }
}
