Use the test data from 10X Genomics `cellranger-arc` program.

The tiny reference and sample FASTQ data can be found in the `external/arc_testrun_files` subfolder of the installation folder (as of v2.0.0).

The sample FASTQs should be placed in the `rawdata` subfolder of the test project folder, and the `essential.vars.groovy` file edited to point to the relevant reference and raw data folders.

Note that the genes GTF file is compressed with cellranger-arc (`genes/genes.gtf.gz`). The `qualimap` tool requires uncompressed GTF files, so an uncompressed version needs to be available for it, and this is the one that should be set in the `essential.vars.groovy` file. The `geneBodyCov2` tool also uses the GTF file, but it can read the compressed version as well.

