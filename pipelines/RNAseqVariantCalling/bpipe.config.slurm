executor="slurm"
queue="short"
commands {
	FastQC {
		walltime="02:00:00"
		procs="4"
		memory="3"
		options="--job-name=FastQC"
	}
	STAR_pe {
      walltime="04:00:00"
		procs="8"
		memory="6"
		options="--job-name=star"
   }
   BAMindexer {
      walltime="01:00:00"
		procs="1"
		memory="1"
		options="--job-name=BAMindexer"
   }
   FilterAndMergeSJtab {
      walltime="01:00:00"
		procs="1"
		memory="1"
		options="--job-name=MergeSJtab"
   }
   GenerateStarIndexFromSJ {
      walltime="04:00:00"
		procs="8"
		memory="32"
		options="--job-name=starIndex"
   }
   STAR_pe_2nd {
      walltime="04:00:00"
		procs="8"
		memory="32"
		options="--job-name=star2nd"
   }
   AddRG {
		walltime="04:00:00"
		procs="1"
		memory="50"
		options="--job-name=AddRG"
	}
   MarkDups {
      walltime="04:00:00"
		procs="1"
		memory="50"
		options="--job-name=MarkDups"
   }
   SplitNCigarReads {
      queue="long"
      procs="8"
		memory="24"
		options="--job-name=SplitNCigarReads"
   }
   BaseRecalibration {
      queue="long"
      procs="8"
		memory="24"
		options="--job-name=BaseRecalibration"
   }
   VariantCallHC {
      queue="long"
      procs="8"
		memory="24"
		options="--job-name=VariantCallHC"
   }
   VariantFiltration {
      walltime="02:00:00"
		procs="8"
		memory="24"
		options="--job-name=VariantFiltration"
   }
}
