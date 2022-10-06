rnaseqc = {
	doc title: "RNA-SeQC: Basic quality control for RNA-seq",
	desc: "efficient RNA-seq quality control and quantification for large cohorts",
	constraints: "",
	author: "Sivarajan Karunanithi"

	output.dir = rnaseqc_vars.outdir
	def RNASEQC_FLAGS = 
		(rnaseqc_vars.legacy ? " --legacy" : "") + 
		(rnaseqc_vars.extra    ? " " + rnaseqc_vars.extra : "")

	def TOOL_ENV = prepare_tool_env("rnaseqc", tools["rnaseqc"]["version"], tools["rnaseqc"]["runenv"])
	def PREAMBLE = get_preamble(stage: stageName, outdir: output.dir, input: new File(input1.prefix).getName())

	//run the chunk
	transform(".bam") to (".bam.gene_reads.gct") {
	exec """
		${TOOL_ENV} &&
		${PREAMBLE} &&
		
		rnaseqc ${rnaseqc_vars.gtf} $input $output.dir ${RNASEQC_FLAGS};
	""","rnaseqc"

	}
}
