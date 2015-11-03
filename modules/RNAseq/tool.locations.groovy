// variables containing the location of the used tools
TOOL_DEPENDENCIES="/fsimb/groups/imb-bioinfocf/common-tools/dependencies/" // your local tools folder
PROJECT_DEPENDENCIES=ESSENTIAL_PROJECT + "/NGSpipe2go/tools/" // please copy the bpipe tools to the project folder and include the location here
TOOL_R=TOOL_DEPENDENCIES + "/R/3.1.2/"
TOOL_FASTQC=TOOL_DEPENDENCIES + "/fastqc/0.11.3"
TOOL_STAR=TOOL_DEPENDENCIES + "/star/2.4.2a/"
TOOL_SAMTOOLS=TOOL_DEPENDENCIES + "/samtools/1.2/samtools"
TOOL_HTSEQ=TOOL_DEPENDENCIES + "/htseq/0.6.1p1"
TOOL_SUBREAD=TOOL_DEPENDENCIES + "/subread/1.4.6-p2/"
TOOL_BEDTOOLS=TOOL_DEPENDENCIES + "/BEDTools/2.16.2/"
TOOL_PICARD=TOOL_DEPENDENCIES + "/picard/1.123/"
TOOL_DUPRADAR=TOOL_DEPENDENCIES + "/imb-forge/dupRadar/0.0.3/"
TOOL_RSeQC=TOOL_DEPENDENCIES + "/RSeQC/2.4/"
TOOL_RNAtypes=TOOL_DEPENDENCIES + "/imb-forge/RNAtypes/"
TOOL_EDGER=PROJECT_DEPENDENCIES + "DE_edgeR/" // this could introduce version differences, if edgeR package is updated. bpipe has no control over the version here.
TOOL_DESeq2=PROJECT_DEPENDENCIES + "DE_DESeq2/" // this could introduce version differences, if DESeq2 package is updated. bpipe has no control over the version here.
TOOL_COLLECT=PROJECT_DEPENDENCIES + "collectBpipeLogs/"
