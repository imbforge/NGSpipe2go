PingPongPro = {
   doc title: "PingPongPro",
         desc: "Calculates the ping-signal for each input transposon (genomic feature)",
         constraints: "Requires PingPongPro http://sourceforge.net/projects/pingpongpro/",
         author: "Antonio Domingues"

   output.dir = PINGPONGPRO_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def FEATURES_NAME = FEATURES_PATH.split("/")[-1].replaceAll(".bed", "")
   def OUTNAME = SAMPLE_NAME + "." + FEATURES_NAME
   def OUT_FOLDER = output.dir + "/" + OUTNAME
   def out1 = OUT_FOLDER + "/ping-pong_signatures.tsv"
   def out2 = OUT_FOLDER + "/transposons.tsv"

   produce(
           out1,
           out2
            ){

      exec """

         module load pingpongpro/${PINGPONGPRO_VERSION} &&

         pingpongpro -i $input -t $FEATURES_PATH -o $OUT_FOLDER

      ""","PingPongPro"
   }
}
