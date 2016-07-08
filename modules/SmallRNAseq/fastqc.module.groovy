FastQC = {
   doc title: "FastQC",
      desc:  "Quality control of input file",
      constraints: "Only supports compressed FASTQ files",
      bpipe_version: "tested with bpipe 0.9.8.7",
      author: "Sergi Sayols; Antonio Domingues"

   output.dir   = FASTQC_OUTDIR
   def FASTQC_FLAGS = "--extract --quiet -t " + FASTQC_THREADS

   transform(".fastq.gz") to ("_fastqc.zip") {
      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(${TOOL_FASTQC} --version | cut -d' ' -f2) 1>&2 &&
         echo '/VERSION INFO' 1>&2 &&

         ${TOOL_FASTQC} $FASTQC_FLAGS -o $output.dir $input
      ""","FastQC"
   }
   forward input
}
