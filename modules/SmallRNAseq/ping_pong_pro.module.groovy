PingPongPro = {
   doc title: "PingPongPro",
         desc: "Calculates the ping-signal for each input transposon (genomic feature)",
         constraints: "Needs PingPongPro in the path http://sourceforge.net/projects/pingpongpro/",
         author: "Antonio Domingues"

   output.dir = PINGPONG_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def OUT_FOLDER = output.dir + "/" + SAMPLE_NAME

   produce(
            OUT_FOLDER + "/ping-pong_signatures.tsv",
            OUT_FOLDER + "/transposons.tsv") {

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         pingpongpro -i $input -t $FEATURES_PATH -o $OUT_FOLDER

      ""","PingPongPro"
   }
}
