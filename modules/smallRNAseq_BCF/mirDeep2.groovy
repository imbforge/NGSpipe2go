miRDeep2 = {
    doc title: "miRDeep2",
        desc:  """Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression. This is step 2""",
        constraints: "Requires mirDeep2.",
        author: "Antonio Domingues, Anke Busch"

    def EXP = input1.split("/")[-1].replaceAll(".arf", "")
    output.dir = miRDeep2_vars.outdir + "/" + EXP

    def TOOL_ENV = prepare_tool_env("mirdeep2", tools["mirdeep2"]["version"], tools["mirdeep2"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".arf", ".fa") to (".tmp") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            reads_fa=`realpath $input2`;
            genome_fa=`realpath $miRDeep2_vars.genome_seq`;
            reads_vs_genome_arf=`realpath $input1`;
            mautre_ref_miRNAs_fa=`realpath $miRDeep2_vars.mature_mirna`;
            mature_other_miRNAs_fa="none";
            hairpin_ref_miRNAs=`realpath $miRDeep2_vars.hairpin_mirna`;

            mkdir -p $output.dir &&
            cd $output.dir &&

            miRDeep2.pl \$reads_fa \$genome_fa \$reads_vs_genome_arf \$mautre_ref_miRNAs_fa \$mature_other_miRNAs_fa \$hairpin_ref_miRNAs -t $miRDeep2_vars.species -c -d -v -r ${EXP} -z ".${EXP}" 2> ${output.dir}/${EXP}.report.log &&
            touch \$(basename $output)
        ""","miRDeep2"
    }
}
