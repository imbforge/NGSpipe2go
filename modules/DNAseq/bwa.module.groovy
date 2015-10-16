//rule for task BWA mem from DNAseq pipeline
//desc: Align single and paired end reads
BWA_pe = {
	doc title: "BWA PE alignment",
	desc:  "Align paired end reads",
	constraints: "Only works with compressed input. Set all global vars.",
	author: "Oliver Drechsel"

	output.dir = MAPPED
    
	// manipulate output file name
	// adapting the solution of Nastasja
	// def OUTPUTFILE = input1
	// OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")
	// produce(OUTPUTFILE + ".bam")
	
	//def SAMPLE_ID = input1
	//SAMPLE_ID = ( SAMPLE_ID =~ /.read_1.fastq.gz/).replaceFirst("") // put suffix to replace in a variable?
	
	//transform('.fastq.gz') to(SAMPLE_ID + '.bam') {
    //transform('.fastq.gz') to('.bam') {
    //transform("bam") {
	transform(".read_1.fastq.gz", ".read_2.fastq.gz") to(".bam") {
	    
        def BWA_FLAGS = "-M -t " + BWA_THREADS //+
                    // " -R @RG\tID:xxx\tSM:xxx\tPL:illumina\tLB:xxx\tPU:xxx" +   + //@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1 
                    // " -p " + BWA_REF
        
        
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
			${TOOL_BWA}/bwa      1>&2 ;
			echo '/VERSION INFO' 1>&2 ;
						
			${TOOL_BWA}/bwa mem $BWA_FLAGS -R \"@RG\\tID:${SAMPLE_NAME}\\tSM:${SAMPLE_NAME}\\tPL:illumina\\tLB:${SAMPLE_NAME}\\tPU:${PLATFORM}\" $ESSENTIAL_BWA_REF $input1 $input2 | ${TOOL_SAMTOOLS} view ${SAMTOOLS_FLAGS} - | ${TOOL_SAMTOOLS} sort -@ ${BWA_THREADS} -O bam -T ${SAMPLE_NAME} -  > ${output} &&
			
			${TOOL_SAMTOOLS} flagstat ${output} 1>&2
			
        ""","BWA_pe"
    }
}