miRDeep2 = {
	doc title: "miRDeep2",
	desc:  """Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression. This is step 2""",
      	constraints: "Requires mirDeep2.",
	author: "Antonio Domingues, Anke Busch"

	def EXP = input1.split("/")[-1].replaceAll(".arf", "")
	output.dir = MIR_OUTDIR + "/" + EXP

	transform(".arf", ".fa") to (".tmp") {

		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_MIRDEEP2}/env.sh &&
 
			mkdir -p $output.dir &&

		        cd $output.dir &&

		        miRDeep2.pl $input2 $GENOME_SEQ $input1 $MATURE_MIRNA none $HAIRPIN_MIRNA -t $SPECIES -c -d -v -r ${EXP} -z ".${EXP}" 2> ${EXP}.report.log &&

			touch ${EXP}.tmp

		""","miRDeep2"
	}
}
