BAMCOVERAGE_CORES="--numberOfProcessors " + Integer.toString(ESSENTIAL_THREADS)
BAMCOVERAGE_OTHER="--smoothLength 60 --binSize 20 --normalizeUsing RPKM --outFileFormat bigwig" //if you want to exclude chromsomes for normalisation e.g. rDNA or mitochondrion add the following parameter --ignoreForNormalization \"chrM, rDNA\", if you like to use offsets, blacklist regions, center reads or anything like it please refer to the deepTools manual, there is even a special modus for Nucleosome detection in Mnase data 
BAMCOVERAGE_OUTDIR=TRACKS
BAMCOVERAGE_FRAGMENTS=ESSENTIAL_FRAGMENT_USAGE
