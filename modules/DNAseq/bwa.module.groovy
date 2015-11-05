//rule for task BWA mem from DNAseq pipeline
//desc: Align single and paired end reads
BWA_pe = {
	doc title: "BWA PE alignment",
	desc:  "Align paired end reads",
	constraints: "Only works with compressed input. Set all global vars.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"

	output.dir = MAPPED
    
	transform(".read_1.fastq.gz", ".read_2.fastq.gz") to(".bam") {
	    
        def BWA_FLAGS = "-M -t " + BWA_THREADS //+
        def SAMTOOLS_FLAGS = "-bhSu"
		
        exec """
        
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_BWA}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
            
			SAMPLE_NAME=\$(basename $output.prefix.prefix) &&

			PLATFORM="genomics" &&
			
			echo 'VERSION INFO'  1>&2 ;
			echo \$(${TOOL_BWA}/bwa 2>&1 | grep Version | cut -d' ' -f2) 1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
						
			${TOOL_BWA}/bwa mem $BWA_FLAGS -R \"@RG\\tID:${SAMPLE_NAME}\\tSM:${SAMPLE_NAME}\\tPL:illumina\\tLB:${SAMPLE_NAME}\\tPU:${PLATFORM}\" $ESSENTIAL_BWA_REF $input1 $input2 | ${TOOL_SAMTOOLS} view ${SAMTOOLS_FLAGS} - | ${TOOL_SAMTOOLS} sort -@ ${BWA_THREADS} -O bam -T ${SAMPLE_NAME} -  > ${output} &&
			
			${TOOL_SAMTOOLS} flagstat ${output} 1>&2
			
        ""","BWA_pe"
    }
}
