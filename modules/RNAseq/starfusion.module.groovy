STAR_Fusion = {
	doc title: "STAR-Fusion",
		desc:  "detection of fusion transcripts from RNA-Seq data",
		constraints: "tab-delimited summary file identifying the fusion pairs",
		bpipe_version: "tested with bpipe 0.9.9",
		author: "Giuseppe Petrosino"

	output.dir = STARFUSION_OUTDIR

    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")

    STARFUSION_FLAGS = STARFUSION_THREADS + " " +
	               STARFUSION_GENOME_LIB

    produce(OUTPUTFILE + "_starfusion.done") {
        exec """
            module load STAR-Fusion/${STARFUSION_VERSION} &&

            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
			
            fi                                       &&

	    STAR-Fusion $STARFUSION_FLAGS --tmpdir $TMPDIR/\$(basename $output.prefix) --left_fq $input1 --right_fq $input2 --output_dir $output.prefix;
            touch $output;		
          
        ""","STAR_Fusion"
    }
}
