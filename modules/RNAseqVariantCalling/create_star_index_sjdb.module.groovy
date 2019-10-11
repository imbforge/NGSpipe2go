// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/create_star_index_sjdb.vars.groovy"

GenerateStarIndexFromSJ = {
    doc title: "GenerateStarIndexFromSJ",
        desc: "Creates a new STAR genome using SJ identified in a previous mapping step. Part of the 2-step mapping. Based on https://code.google.com/p/rna-star/issues/detail?id=7",
        constraints: "STAR STAR_2.4.2a",
        author: "Antonio Domingues"

    output.dir = GenerateStarIndexFromSJ_vars.outdir

    def STAR_FLAGS =
        " --runMode genomeGenerate " +
        " --readFilesCommand zcat "  +
        (GenerateStarIndexFromSJ_vars.outdir_2nd_index ? " --genomeDir "              + GenerateStarIndexFromSJ_vars.outdir_2nd_index : "") +
        (GenerateStarIndexFromSJ_vars.genome_ref       ? " --genomeFastaFiles "       + GenerateStarIndexFromSJ_vars.genome_ref : "") +
        (GenerateStarIndexFromSJ_vars.sjdbfile         ? " --sjdbFileChrStartEnd "    + GenerateStarIndexFromSJ_vars.sjdbfile   : "") +
        (GenerateStarIndexFromSJ_vars.overhang         ? " --sjdbOverhang "           + GenerateStarIndexFromSJ_vars.overhang   : "") +
        (GenerateStarIndexFromSJ_vars.threads          ? " --runThreadN "             + GenerateStarIndexFromSJ_vars.threads    : "") +
        (GenerateStarIndexFromSJ_vars.maxram           ? " --limitGenomeGenerateRAM " + GenerateStarIndexFromSJ_vars.maxram     : "") +
        (GenerateStarIndexFromSJ_vars.bufsize          ? " --limitIObufferSize "      + GenerateStarIndexFromSJ_vars.bufsize    : "") +
        (GenerateStarIndexFromSJ_vars.extra            ? " "                          + GenerateStarIndexFromSJ_vars.extra      : "")

    def TOOL_ENV = prepare_tool_env("star", tools["star"]["version"], tools["star"]["runenv"])
    def PREAMBLE = get_preamble("GenerateStarIndexFromSJ")

    produce("sjdbInfo.txt"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            STAR $STAR_FLAGS
        ""","GenerateStarIndexFromSJ"
    }
}
