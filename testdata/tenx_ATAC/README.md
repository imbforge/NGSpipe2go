Use the test data from 10X Genomics `cellranger-atac` program.

The tiny reference and sample FASTQ data can be found in the `external/atac_testrun_ref` and `external/cellranger_atac_tiny_fastq` subfolders of the installation folder (as of v2.0.0).

The sample FASTQs should be placed in the `rawdata` subfolder of the test project folder, and the `essential.vars.groovy` file edited to point to the relevant reference and raw data folders.

Note that the genes GTF file is compressed with cellranger-atac (`genes/genes.gtf.gz`). This is not important as no NGSpipe2go pipeline modules use it in this pipeline.

