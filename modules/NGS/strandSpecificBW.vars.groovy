//VARS for task strandspecific
STRANDSPECIFICBIGWIG_CORES=4
STRANDSPECIFICBIGWIG_OTHER="--smoothLength 60 --binSize 20 --normalizeUsingRPKM --skipNonCoveredRegions --outFileFormat bedgraph" // if you use > v3 modify the vars to --normalizeUsing RPKM since the API of deeptools changed
STRANDSPECIFICBIGWIG_OUTDIR=TRACKS + "/strandspecific"
