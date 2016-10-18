// variables containing the location of the used tools
TOOL_DEPENDENCIES="/opt/" // your local tools folder
PROJECT_DEPENDENCIES= ESSENTIAL_PROJECT + "/NGSpipe2go/tools/" // please copy the NGSpipe2Go tools to the project folder and include the location here
TOOL_JAVA=TOOL_DEPENDENCIES + "/jdk/1.8.0_102/jre/bin/java"
TOOL_SAMTOOLS=TOOL_DEPENDENCIES + "/samtools/1.3/samtools"
TOOL_FASTQC=TOOL_DEPENDENCIES + "/fastqc/0.11.3"
TOOL_BWA=TOOL_DEPENDENCIES + "/bwa/0.7.12"
TOOL_PICARD=TOOL_DEPENDENCIES + "/picard/2.7.0/picard.jar"
TOOL_GATK=TOOL_DEPENDENCIES + "/GATK/GATK-3.4-46"
TOOL_COLLECT=PROJECT_DEPENDENCIES + "/collectBpipeLogs/" // this is provided by NGSpipe2go and hence locates to the project folder
