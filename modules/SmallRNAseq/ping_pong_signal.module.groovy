load MODULE_FOLDER + "SmallRNAseq/ping_pong_signal.vars.groovy"

PingPongSignal = {
   doc title: "Ping-Pong signal",
         desc:  "Calculates the 10bp overlap frequency of read pairs, known as ping-pong. Outputs a plot with the signal and z-score.",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = PINGPONG_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def FEATURES_NAME = FEATURES_PATH.split("/")[-1].replaceAll(".bed", "")
   def OUTNAME = SAMPLE_NAME + "." + FEATURES_NAME
   def OUT_FOLDER = output.dir + "/" + OUTNAME
   def OUT_PDF = OUT_FOLDER + "/figure/" + OUTNAME + ".ppPlot.pdf"

   produce(
      OUT_PDF
      ) {

      exec """
         module load R/${R_VERSION} &&

         python $PINGPONG_TOOL_PATH -b $input --minsize 19 --maxsize 35 --primary sense --outFolder $OUT_FOLDER --intervals $FEATURES_PATH

      ""","PingPongSignal"
   }
}
