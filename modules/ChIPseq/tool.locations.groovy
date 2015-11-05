// variables containing the location of the used tools
TOOL_DEPENDENCIES="/opt/" // your local tools folder
PROJECT_DEPENDENCIES=ESSENTIAL_PROJECT + "/NGSpipe2go/tools/" // please copy the NGSpipe2Go tools to the project folder and include the location here
TOOL_R=TOOL_DEPENDENCIES + "/R/3.2.2/"
TOOL_FASTQC=TOOL_DEPENDENCIES + "/fastqc/0.11.3"
TOOL_BOWTIE=TOOL_DEPENDENCIES + "bowtie/1.1.1/"
TOOL_SAMTOOLS=TOOL_DEPENDENCIES + "/samtools/1.2/bin/samtools"
TOOL_BEDTOOLS=TOOL_DEPENDENCIES + "/BEDTools/2.16.2/"
TOOL_PICARD=TOOL_DEPENDENCIES + "/picard/1.123/"
TOOL_RSeQC=TOOL_DEPENDENCIES + "/RSeQC/2.4/"
TOOL_ENCODEqc=PROJECT_DEPENDENCIES + "/encodeChIPqc/"
TOOL_MACS2=TOOL_DEPENDENCIES + "/macs2/2.1.0/"
TOOL_COLLECT=PROJECT_DEPENDENCIES + "/collectBpipeLogs/"
