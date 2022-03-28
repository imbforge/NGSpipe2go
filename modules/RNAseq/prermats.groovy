PRERMATS = {
    doc title: "rMATS preprocessor for the target and contrast file",
    desc: "rMATS Preprocessing of targets and contrast file",
    constraints: "",
    author: "Nastasja Kreim, Modified by: Sivarajan Karunanithi"

    output.dir = PRERMATS_vars.outdir
    //read in the contrasts file 
    def contrasts = new ArrayList()
    def f = new File(PRERMATS_vars.contrasts)
    f.eachLine { line, number -> 
        if(number == 1) { return } // skip headerline of contrasts file
        contrasts.add(line) 
    }
    //extract only the contrasts name and 
    //create respective targets file name for rmats
    contrasts.eachWithIndex { item, index ->
        def item_list = item.split("\t")
        contrasts[index]=(item_list[1]).substring(1,item_list[1].length()-1) + PRERMATS_vars.suffix
        println contrasts[index]

        // In case of pre-existing results, the time stamps of all inputs have
        // to be compared to those of the outputs. While this should be done by 
        // bpipe automatically, it does not work correctly in case of a
        // "produce(array)" statement to generate the output(s). In case of 
        // pre-existing outputs, new outputs are written, but their time stamp is
        // not updated fast enough before the next module is started. Thus, the 
        // next module will see the old time stamp and, thus, not re-create its 
        // own ouputs. In order to correct for this, the time stamps are checked 
        // manually and possibly pre-existing outputs are deleted:
        File f_output = new File(output.dir + "/" + contrasts[index])
        if(f_output.exists()) {
            def newcount=0
            inputs.each { val -> 
                File f_input = new File(val.toString())
                if(f_input.lastModified() > f_output.lastModified()) { newcount++ }
            }
            if(newcount > 0) {
                println contrasts[index] + " was already there and deleted " + (f_output.delete() ? "successfully" : "unsuccessfully")
            }
        }
    }

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    // Todo: handle case when multiple contrast name exists for same contrasts. 
    // Eg. This may happen if one explores different mmatrix values for same contrasts.
    produce(contrasts) {
        exec """

            ${PREAMBLE} &&

            for i in `cut -f2 $PRERMATS_vars.contrasts | sort | uniq`; do 
                if [[ "\${i}" == "contrast" ]];
                then
                   echo "Skipping header line of contrasts file\n";
                   continue;
                fi;
                contrast=\${i:1:-1};
                groups[0]=`echo \${i:1:-1} | cut -d- -f1`;
                groups[1]=`echo \${i:1:-1} | cut -d- -f2`;
                groups=`echo \${groups[0]} \${groups[1]}`;
                echo \${groups};
                awk -v G="\$groups" 'BEGIN {OFS="\t"}{ split(G, g, / /)} { if( \$3 == g[1] || \$3 == g[2] ) { split(\$2,f,"."); \$2=f[1]; print \$0 } }' $PRERMATS_vars.targets > $output.dir/\${contrast}_targets_rMATS.txt;
            done
        """, "PRERMATS"
    }
}
