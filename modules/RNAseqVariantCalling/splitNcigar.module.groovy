// rule for task GATK SplitNCigarReads from the RNAseq variant calling pipeline
SplitNCigarReads = {
   doc title: "GATK SplitNCigarReads",
       desc: "Splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions",
       constraints: "GATK version >= 3.5",
       author: "Antonio Domingues"

   output.dir = OUTDIR_STAR2ND

   def JAVA_FLAGS = "-Xmx" + SPLITNCIGAR_MAXMEM
   def SPLITCIGAR_FLAGS  = " -R " + GATK_REF +
                         " -rf " + READ_FILTER_FLAG +
                         " -RMQF " + MAP_Q_FROM_FLAG +
                         " -RMQT " + MAP_Q_TO_FLAG +
                         " -U " + UNSAFE_FLAG

   transform (".duprm.bam") to (".duprm.split.bam"){

      exec """
            echo 'VERSION INFO'  1>&2 &&
            echo \$(java -jar ${TOOL_GATK}/GenomeAnalysisTK.jar --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            java $JAVA_FLAGS -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T SplitNCigarReads -I $input -o $output $SPLITCIGAR_FLAGS
      ""","SplitNCigarReads"
   }
}
