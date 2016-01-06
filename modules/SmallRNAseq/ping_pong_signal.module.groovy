PingPongSignal = {
   doc title: "Ping-Pong signal",
         desc:  "Calculates the 10bp overlap frequency of read pairs, known as ping-pong. Outputs a plot with the signal and z-score.",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = PINGPONG_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def OUT_FOLDER = output.dir + "/" + SAMPLE_NAME

   produce(OUT_FOLDER + "/figure/pp_freq.txtppPlot.pdf") {

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&
         source ${TOOL_PYTHONENV} &&

         python $PINGPONG_TOOL_PATH -b $input --minsize 19 --maxsize 35 --primary sense --outFolder $OUT_FOLDER --intervals $FEATURES_PATH

      ""","PingPongSignal"
   }
}
