pbc = {
    doc title: "PBC",
        desc:  "PCR Bottleneck Coefficient",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = QC + "/pbc"

    transform(".bam") to("_PBC.csv") {
        exec """
            module load R/${R_VERSION} &&

            Rscript ${TOOL_ENCODEqc}/PBC.R $input && mv ${input.prefix}_PBC.csv $output.dir
        ""","pbc"
    }

    forward input

}
