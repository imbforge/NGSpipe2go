// rule to add read groups to bam
@filter('rg')
AddReadGroup = {
   dic title: "AddReadGroup",
   desc: "Adds reads groups to bam as part of the GATK pipeline",
   constraints: "requires picard tools"
   author: "Antonio Domingues"

   output.dir = OUTDIR_STAR2ND
   DEF JAVA_FLAGS = "-Xmx" + RG_MAXMEM

   exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         SAMPLE_NAME=\$(basename $input.prefix) &&
         PLATFORM="genomics" &&

         ${TOOL_JAVA}/java ${JAVA_FLAGS} -jar ${TOOL_PICARD} AddOrReplaceReadGroups I=$input O=$output SO=coordinate RGID=${SAMPLE_NAME} RGLB=${SAMPLE_NAME} RGPL=illumina RGPU=${PLATFORM} RGSM=${SAMPLE_NAME}

   ""","AddRG"
}
