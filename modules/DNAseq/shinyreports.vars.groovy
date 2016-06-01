//vars for task shinyReports from catalog miscellaneous, version 1
SHINYREPS_PROJECT=PROJECT	//project directory
SHINYREPS_LOG=LOGS			//where the logs lie
SHINYREPS_QC=QC				//where the QC lie
SHINYREPS_RES=RESULTS		//where the results lie
SHINYREPS_PREFIX=ESSENTIAL_SAMPLE_PREFIX	//standard sample prefix
SHINYREPS_FASTQC_OUT=FASTQC_OUTDIR		//where the Fastqc output lie
SHINYREPS_FASTQC_LOG=LOGS + "/FastQC"		//where the Fastqc logs lie
SHINYREPS_BWA_LOG=LOGS + "/BWA_pe"	//where the BWA (samtools flagstat) logs lie
SHINYREPS_BWA_SUFFIX=".bam.log"	//extension given to the BWA log files
SHINYREPS_MARKDUPS_LOG=LOGS + "/MarkDups"	//where the picard MarkDuplicates logs lie
SHINYREPS_GATKug_LOG=LOGS + "/VariantCallUG"	//where the GATK UnifiedGenotyper logs lie
SHINYREPS_GATKug_SUFFIX=".UG.vcf.gz.log"	//extension given to the GATK UnifiedGenotyper log files
SHINYREPS_GATKhc_LOG=LOGS + "/VariantCallHC"	//where the GATK HaplotypeCaller logs lie
SHINYREPS_GATKhc_SUFFIX=".HC.vcf.gz.log"	//extension given to the GATK Haplotypecaller log files
SHINYREPS_GATKvarianteval=QC + "/GATK_varianteval" // location of GATK variantEval results
SHINYREPS_GATKvarianteval_SUFFIX=".report"	//extension given to the GATK VariantEval output files
SHINYREPS_RES_GATKhc_SUFFIX=".HC.vcf.gz" // extension of the GATK HaplotypeCaller output files

