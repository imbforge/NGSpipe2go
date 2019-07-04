//rule for task merge_bams
//desc: Counting reads in features with featureCounts
MergeBam = {
   doc title: "Merge bam files",
   desc:  "Merges bam files following any given pipeline defined pattern",
   constraints: "Unless modified, the name for the merged replicates will be determined by removing the pattern _rep[1-9] from the name of the first input. Change the code bellow if the pattern of your samples is different",
   bpipe_version: "tested with bpipe 0.9.9.3",
   author: "Antonio Domingues"

   output.dir  = MERGEDBAMS_OUTDIR
   def EXP = input1.split("/")[-1].replaceAll(".bam", "").replaceAll("_rep\\d+", "")

   // run the chunk
   produce(EXP + ".merged.bam") {
      exec """
         module load samtools/${SAMTOOLS_VERSION} &&

         echo $inputs &&
         samtools merge $output $inputs 

      ""","MergeBam"
   }
}
