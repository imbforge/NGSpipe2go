// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mirDeep2.vars.groovy"

miRDeep2 = {
    doc title: "miRDeep2",
        desc:  """Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression. This is step 2""",
        constraints: "Requires mirDeep2.",
        author: "Antonio Domingues, Anke Busch"

    def EXP = input1.split("/")[-1].replaceAll(".arf", "")
    output.dir = MIR_OUTDIR + "/" + EXP

    def TOOL_ENV = prepare_tool_env("mirdeep2", tools["mirdeep2"]["version"], tools["mirdeep2"]["runenv"])

    transform(".arf", ".fa") to (".tmp") {
        exec """
            ${TOOL_ENV} &&

            mkdir -p $output.dir &&
            cd $output.dir &&
            miRDeep2.pl $input2 $GENOME_SEQ $input1 $MATURE_MIRNA none $HAIRPIN_MIRNA -t $SPECIES -c -d -v -r ${EXP} -z ".${EXP}" 2> ${EXP}.report.log &&
            touch ${EXP}.tmp
        ""","miRDeep2"
    }
}
