PingPongSignal = {
   doc title: "Ping-Pong signal",
         desc:  "Calculates the 10bp overlap frequency of read pairs, known as ping-pong. Outputs a plot with the signal and z-score.",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = PINGPONG_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def OUT_FOLDER = output.dir + "/" + SAMPLE_NAME
   "figure/Tdrd9-het-Trunk-zili-W3.ppPlot.pdf"
   produce(OUT_FOLDER + "/figure/" + SAMPLE_NAME + ".ppPlot.pdf") {

      exec """
         source ${TOOL_PYTHONENV} &&

         python $PINGPONG_TOOL_PATH -b $input --minsize 19 --maxsize 35 --primary sense --outFolder $OUT_FOLDER --intervals $FEATURES_PATH

      ""","PingPongSignal"
   }
}
