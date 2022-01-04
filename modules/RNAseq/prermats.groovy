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
        item = ( item =~ /\t.*$/).replaceFirst("")
        contrasts[index]=item + PRERMATS_vars.suffix
        println contrasts[index]
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
                contrast=`grep \${i} $PRERMATS_vars.contrasts | cut -f1`;
                groups[0]=`echo \${i:1:-1} | cut -d- -f1`;
                groups[1]=`echo \${i:1:-1} | cut -d- -f2`;
                groups=`echo \${groups[0]} \${groups[1]}`;
                echo \${groups[0]} \${groups[1]};
                awk -v G="\$groups" 'BEGIN { split(G, g, / /)} { if( \$3 == g[1] || \$3 == g[2] ) print \$0 }' $PRERMATS_vars.targets > $output.dir/\${contrast}_targets_rMATS.txt;
            done
        """, "PRERMATS"
    }
}
