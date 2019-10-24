// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/star.vars.groovy"

STAR = {
    doc title: "STAR alignment",
        desc:  "Align single/paired end reads",
        constraints: "Paired end reads expected to have .R1 and .R2 suffix.",
        bpipe_version: "tested with bpipe 0.9.9",
        author: "Sergi Sayols"

    output.dir = STAR_vars.outdir

    // create the LOGS/STAR folder if it doesn't exists
    def F_LOG = new File(STAR_vars.logdir)
    if(! F_LOG.exists()) {
        F_LOG.mkdirs()
    }

    // calculate name of the sample being processed (added paired end support)
    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    // star flags
    def STAR_FLAGS =
        " --runMode alignReads "        +
        " --genomeLoad NoSharedMemory " +
        " --outStd SAM "                +
        " --outSAMattributes Standard " +
        " --outSJfilterReads Unique "   +
        " --readFilesCommand zcat "     +
        " --outFileNamePrefix "         + STAR_vars.logdir + "/" + OUTPUTFILE +
        (STAR_vars.unmapped_bam ? " --outSAMunmapped " + STAR_vars.unmapped_bam   : "") +
        (STAR_vars.ref          ? " --genomeDir "      + STAR_vars.ref            : "") +
        (STAR_vars.threads      ? " --runThreadN "     + STAR_vars.threads        : "") +
        (STAR_vars.mm           ? " --outFilterMismatchNoverLmax " + STAR_vars.mm : "") +
        (STAR_vars.overhang     ? " --sjdbOverhang "   + STAR_vars.overhang       : "") +
        (STAR_vars.gtf          ? " --sjdbGTFfile "    + STAR_vars.gtf            : "") +
        (STAR_vars.extra        ? " "                  + STAR_vars.extra          : "")

    // samtools flags
    def SAMTOOLS_VIEW_FLAGS = "-bhSu" +
        (STAR_vars.filter_sec ? " -F 256" : "")
    def SAMTOOLS_SORT_FLAGS = " -O bam" +
        (STAR_vars.samtools_threads ? " -@ " + STAR_vars.samtools_threads : "")

    def TOOL_ENV = prepare_tool_env("star", tools["star"]["version"], tools["star"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("STAR")

    // code chunk
    // TODO: warn if the genome index was created using another version of STAR?
    produce(OUTPUTFILE + ".bam", OUTPUTFILE + "Log.final.out") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            if [ -e \${TMP}/$OUTPUTFILE ]; then
                echo 'removing old STAR tmp folder';
                rm -r \${TMP}/$OUTPUTFILE*;
            fi &&

            STAR $STAR_FLAGS --outTmpDir \${TMP}/$OUTPUTFILE --readFilesIn $inputs | samtools view $SAMTOOLS_VIEW_FLAGS - | samtools sort $SAMTOOLS_SORT_FLAGS -T \${TMP}/${OUTPUTFILE}_sort - > $output1 &&

            mv ${STAR_vars.logdir}/${OUTPUTFILE}SJ.out.tab $output.dir &&
            ln -s ${STAR_vars.logdir}/${OUTPUTFILE}Log.final.out $output.dir
        ""","STAR"
    }
}
