// rule to add read groups to bam
AddRG = {
   doc title: "AddReadGroup",
   desc: "Adds reads groups to bam as part of the GATK pipeline",
   constraints: "Picard tools version >= 1.141"
   author: "Antonio Domingues"

   output.dir = OUTDIR_STAR2ND
   def JAVA_FLAGS = "-Xmx" + RG_MAXMEM
   def EXP = input1.split("/")[-1].replaceAll(".bam", "")

   transform(".bam") to (".rg.bam"){
      exec """
            echo 'VERSION INFO'  1>&2 &&
            echo \$(java -jar ${TOOL_PICARD}/picard.jar AddOrReplaceReadGroups --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            PLATFORM="genomics" &&

            java $JAVA_FLAGS -jar ${TOOL_PICARD}/picard.jar AddOrReplaceReadGroups I=$input O=$output SO=coordinate RGID=${EXP} RGLB=${EXP} RGPL=illumina RGPU=${PLATFORM} RGSM=${EXP}

      ""","AddRG"
   }
}
