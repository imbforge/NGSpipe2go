PingPongPro = {
   doc title: "PingPongPro",
         desc: "Calculates the ping-signal for each input transposon (genomic feature)",
         constraints: "Requires PingPongPro http://sourceforge.net/projects/pingpongpro/",
         author: "Antonio Domingues"

   output.dir = PINGPONG_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def OUT_FOLDER = output.dir + "/" + SAMPLE_NAME

   produce(
            OUT_FOLDER + "/ping-pong_signatures.tsv",
            OUT_FOLDER + "/transposons.tsv") {

      exec """

         module load pingpongpro/${PINGPONGPRO_VERSION} &&

         pingpongpro -i $input -t $FEATURES_PATH -o $OUT_FOLDER

      ""","PingPongPro"
   }
}
