assignSouporcellCluster = {
    doc title: "Assign souporcell clusters after running the demux_gt module",
        desc:  "Determine which clusters are the shared clusters between souporcell runs with shared samples (sub-samples)",
        constraints: "This does not yet determine which cluster refers to which sub-dample in the sample sheet. For this we would need externally provided genotypes of the subsamples.",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = assignSouporcellCluster_vars.outdir + "/" 


    def TOOL_ENV = prepare_tool_env("souporcell", tools["souporcell"]["version"], tools["souporcell"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(output.dir + "/assignSouporcellCluster.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            headerposFile=\$(head -n1 ${assignSouporcellCluster_vars.targets} | tr "\\t" "\\n" | grep -nx file | cut -d":" -f1) &&
            filenames=(\$(tail -n +2 ${assignSouporcellCluster_vars.targets} | cut -f\${headerposFile} | sort | uniq | sed 's/.fastq.gz*\$//' | sed 's/_S[0-9]*_L[0-9]*_R[12]*_00[1-4]*//')) &&
            totalFiles=\${#filenames[*]} &&

            for ((i=0; i<\${totalFiles}-1; i++)); do
                for ((j=\$i+1; j<\$totalFiles; j++)); do 
                    echo \${i}_\${j};
                    echo \${filenames[\$i]}_vs_\${filenames[\$j]};
                    sample1=${assignSouporcellCluster_vars.souporcelldir}/\${filenames[\$i]};
                    sample2=${assignSouporcellCluster_vars.souporcelldir}/\${filenames[\$j]};
                    echo \$sample1;
                    echo \$sample2;

                    clusterS1=\$(cat ${assignSouporcellCluster_vars.targets} | cut -f\${headerposFile} | grep \${filenames[\$i]} | wc -l);
                    clusterS2=\$(cat ${assignSouporcellCluster_vars.targets} | cut -f\${headerposFile} | grep \${filenames[\$j]} | wc -l);
                    sharedCluster=\$(( \$clusterS1 < \$clusterS2 ? \$clusterS1 : \$clusterS2 ));
                    echo sharedCluster_\$sharedCluster;

                    shared_samples.py -1 \$sample1 -2 \$sample2 -n \$sharedCluster > ${output.dir}/\${filenames[\$i]}_vs_\${filenames[\$j]}.txt;
                done 
            done &&

            touch $output

        ""","assignSouporcellCluster"
    }
}



