Use the test data from 10X Genomics `cellranger` program.

The tiny reference and sample FASTQ data can be found in the `external/cellranger_tiny_ref` and `external/cellranger_tiny_fastq` subfolders of the installation folder (as of v6.0.0).

The sample FASTQs should be placed in the `rawdata` subfolder of the test project folder, and the `essential.vars.groovy` file edited to point to the relevant reference and raw data folders.

Note that for the test data associated with cellranger (v6.0.0), the `ESSENTIAL_FEATURETYPE` variable should be set to the (ENSEMBL-associated) `gene_biotype`, not the (GENCODE-associated) `gene_type`. The full human and mouse 10X datasets, as well as the test datasets for cellranger-atac and cellranger-arc, use the GENCODE-style `gene_type`.

