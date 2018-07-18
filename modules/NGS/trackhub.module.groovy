trackhub = {
    doc title: "Trackhub",
        desc:  "Generate UCSC track hub to display project tracks",
        constraints: "Uses configuration file, which should have been generated earlier",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Martin Oti"

    def trackhub_FLAGS = "ESSENTIAL_PROJECT=" + ESSENTIAL_PROJECT + " " +
        "TRACKHUB_CONFIG=" + TRACKHUB_CONFIG

    transform(".yaml") to (".done") {
        exec """
           module load kentUtils/${KENTUTILS_VERSION} &&
           module load R/${R_VERSION} &&
           touch $output;
           Rscript ${TOOL_TRACKHUB}/Make_Trackhub.R $trackhub_FLAGS &&
           if [ \$? -ne 0 ]; then rm $output; fi;
        ""","trackhub"
    }
}

