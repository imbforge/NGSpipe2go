//rule for task STAR_se from catalog RNAseq, version 1
//desc: Align single end reads
STAR_se = {
	doc title: "STAR SE alignment",
		desc:  "Align single end reads",
		constraints: "Only works with compressed input. Set all global vars.",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = MAPPED

	// create the TMP folder if it doesn't exists
	F_TMP = new File(TMP)
	if(! F_TMP.exists()) { 
		F_TMP.mkdirs()
	}
	// create the LOGS/STAR folder if it doesn't exists
	F_LOG = new File(LOGS + "/STAR_se")
	if(! F_LOG.exists()) {
		F_LOG.mkdirs()
	}
	
	def OUTPUTFILE = input1
	int path_index = OUTPUTFILE.lastIndexOf("/")
	OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
	println(OUTPUTFILE)
        if( ESSENTIAL_PAIRED == "yes") {
		OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")
	}
	else {
		OUTPUTFILE = (OUTPUTFILE =~ /.fastq.gz/).replaceFirst("")
	}

	// code chunk
	produce(OUTPUTFILE + ".bam", OUTPUTFILE + "Log.final.out") {
//	transform("_R1.fastq.gz") to(".bam", "Log.final.out") {
		// flags
		def SAMPLE = new File(input.prefix.prefix)
		def STAR_FLAGS = "--runMode alignReads "        +
                         "--genomeLoad NoSharedMemory " +
                         "--outStd SAM "                +
                         "--outSAMattributes Standard " +
                         "--outSJfilterReads Unique "   +
                         "--readFilesCommand zcat "     +
                         "--outFileNamePrefix " + LOGS + "/STAR_se/" + SAMPLE.name + " " +
                         "--outTmpDir " + TMP + "/" + SAMPLE.name + " " +
						 STAR_UNMAPPED_BAM + " " +
						 STAR_UNMAPPED_OUT + " " +
                         STAR_MAXRAM   + " " +
                         STAR_BUFSIZE  + " " +
                         STAR_REF      + " " +
                         STAR_THREADS  + " " +
                         STAR_MM       + " " +
                         STAR_MULTIMAP + " " +
                         STAR_MININTRO + " " +
                         STAR_OVERHANG + " " +
                         STAR_GTF      + " " +
                         STAR_EXTRA

		def SAMTOOLS_FLAGS = "-bhSu"
		if(STAR_FILTER_SEC == "YES") {
			SAMTOOLS_FLAGS = " -F 256 " + SAMTOOLS_FLAGS
		}
		
		def MMR_FLAGS = "-F 0 " + // --
				"-b " +
				"-B " +
				MMR_THREADS + " " +
				MMR_WINDOWSIZE + " " +
				MMR_ANNOTATION + " " +
				MMR_EXTRA
		
		
		// TODO: change to latest or at least try to warn, if the genome index was created using the wrong version of STAR
		// DONE: replace ""source ${TOOL_DEPENDENCIES}/star/default/env.sh""
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_STAR}/env.sh &&
			source ${TOOL_MMR}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			else 
				export TMPDIR=/tmp;
			fi                                          &&
			
			if [ -e $TMP/$SAMPLE.name ];
			then
				echo 'removing old STAR tmp folder';
				rm -r $TMP/$SAMPLE.name*;
			fi &&
			
			echo 'VERSION INFO'  1>&2 &&
			echo \$(STAR --version) 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&
			
			STAROUTPUTFILE=\$(basename $output1.prefix) &&
			
			STAR $STAR_FLAGS --readFilesIn $inputs | ${TOOL_SAMTOOLS} view $SAMTOOLS_FLAGS - | ${TOOL_SAMTOOLS} sort -n $STAR_SAMTOOLS_THREADS - > \${TMPDIR}/\${STAROUTPUTFILE}_beforeMMR.bam &&
			echo "STAR done." &&
			
			mv ${LOGS}/STAR_se/${SAMPLE.name}SJ.out.tab $output.dir &&
			ln -s ${LOGS}/STAR_se/${SAMPLE.name}Log.final.out $output.dir &&
			
			mmr $MMR_FLAGS -o \${TMPDIR}/\${STAROUTPUTFILE}_afterMMR.bam \${TMPDIR}/\${STAROUTPUTFILE}_beforeMMR.bam && 
			echo "MMR done." &&
			
			${TOOL_SAMTOOLS} view -b -f 4 \${TMPDIR}/\${STAROUTPUTFILE}_beforeMMR.bam | ${TOOL_SAMTOOLS} merge -np \${TMPDIR}/\${STAROUTPUTFILE}_afterMMRplusUnmapped.bam \${TMPDIR}/\${STAROUTPUTFILE}_afterMMR.bam - &&

			${TOOL_SAMTOOLS} sort \${TMPDIR}/\${STAROUTPUTFILE}_afterMMRplusUnmapped.bam > $output1
		""","STAR_se"
	}
}
