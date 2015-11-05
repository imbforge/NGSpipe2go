//rule for task GATK realignment from DNAseq pipeline
//desc: Indel realignment provided by GATK
IndelRealignment = {
    
    doc title: "GATK IndelRealignment",
	desc:  "Realign BAM files at Indel positions, using GATK",
	constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
	author: "Oliver Drechsel"

	output.dir = MAPPED
    
    transform (".bam") to (".realignment.targets.bed", ".realigned.bam") {
    // filter("realigned") {
        
        def GATK_FLAGS = "-known gold_indels.vcf "
        
        // check if a region limit was provided
        if (ESSENTIAL_CALL_REGION!=null && ESSENTIAL_CALL_REGION.length()>0) {
            
            GATK_FLAGS = GATK_FLAGS + " -L " + ESSENTIAL_CALL_REGION
            
        } else {
            
            GATK_FLAGS = ""
            
        }
        
        exec """
            export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi                                          &&
        
            java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $GATK_THREADS -R $ESSENTIAL_BWA_REF -I $input -o $output1 &&
        
            java -Djava.io.tmpdir=$TMPDIR -jar $TOOL_GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $ESSENTIAL_BWA_REF -I $input -targetIntervals $output1 -o $output2 
        
        ""","IndelRealignment"
    }
    
}
