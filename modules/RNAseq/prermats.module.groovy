load MODULE_FOLDER + "RNAseq/prermats.vars.groovy"

PREMATS = {
    doc title: "rMATs prepocessor for the target and contrast file",
    desc: "rMATS Prepocessing of targets and contrast file",
    constraints: "",
    author: "Nastasja Kreim"

    output.dir = PREMATS_OUTDIR
    //read in the contrasts file 
    def contrasts = new ArrayList()
    def f = new File(PREMATS_CONTRASTS)
    f.eachLine { line -> contrasts.add(line) }
    //try to replace the contrast lines by the first part
    contrasts.eachWithIndex { item, index -> 
        //we have to skipt the first line
        item = ( item =~ /=.*$/).replaceFirst("")
            println item
            contrasts[index]=item + PREMATS_POSTFIX
            println contrasts[index]
    }

    produce(contrasts) {
        exec """
            for i in `cat $PREMATS_CONTRASTS`; do
                groups=(\${i//=/" "});
                contrast=\${groups[0]};
                tmp=\${groups[1]};
                groups=(\${tmp//-/ });
                groups[0]=`echo \${groups[0]} | sed 's/(//g'`;
                groups[1]=`echo \${groups[1]} | sed 's/)//g'`;
                groups=`echo \${groups[0]} \${groups[1]}`;
                echo \${groups[0]} \${groups[1]};
                awk -v G="\$groups" 'BEGIN { split(G, g, / /)} { if( \$2 == g[1] || \$2 == g[2] ) print \$0 }' $PREMATS_TARGETS > $output.dir/\${contrast}_targets_rMats.txt;
            done
        """, "PREMATS"
    }
}
