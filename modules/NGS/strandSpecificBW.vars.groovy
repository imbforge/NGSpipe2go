//VARS for task strandspecific
STRANDSPECIFICBIGWIG_CORES=4
STRANDSPECIFICBIGWIG_OTHER="--smoothLength 60 --binSize 20 --normalizeUsingRPKM --skipNonCoveredRegions --outFileFormat bedgraph"
STRANDSPECIFICBIGWIG_OUTDIR=TRACKS + "/strandspecific"
STRANDSPECIFICBIGWIG_CHRSIZES="/fsimb/groups/imb-bioinfocf/projects/roignant2/imb_roignant_2016_meta_02_JA_ChIP_4SU_RNA/reference/translation_table_dm6/dm6.chrom.sizes"
STRANDSPECIFICBIGWIG_TRANSTABLE="/fsimb/groups/imb-bioinfocf/projects/roignant/imb_roignant_2016_02_ND_RNA/reference/translation_table_dm6/dm6_translation_table.txt"
