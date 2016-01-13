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
            if [ -n "\$LSB_JOBID" ]; then
               export TMPDIR=/jobdir/\${LSB_JOBID};
            fi &&

            echo 'VERSION INFO'  1>&2 &&
            echo \$(${TOOL_JAVA}/java -jar ${TOOL_PICARD} AddOrReplaceReadGroups --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            PLATFORM="genomics" &&

            ${TOOL_JAVA}/java ${JAVA_FLAGS} -jar ${TOOL_PICARD} AddOrReplaceReadGroups I=$input O=$output SO=coordinate RGID=${EXP} RGLB=${EXP} RGPL=illumina RGPU=${PLATFORM} RGSM=${EXP}

      ""","AddRG"
   }
}
