// variables containing the location of the used tools
TOOL_DEPENDENCIES="/fsimb/groups/imb-kettinggr/common_bin" // your local tools folder
PROJECT_DEPENDENCIES= ESSENTIAL_PROJECT + "/NGSpipe2go/tools/" // please copy the NGSpipe2Go tools to the project folder and include the location here
TOOL_JAVA="/usr/lib/jvm/java-7-openjdk-amd64/jre/bin/"
TOOL_SAMTOOLS=TOOL_DEPENDENCIES + "/samtools/1.2/bin/samtools"
TOOL_FASTQC=TOOL_DEPENDENCIES + "/FastQC/0.11.4/fastqc"
TOOL_BWA=TOOL_DEPENDENCIES + "/bwa/0.7.12"
TOOL_PICARD=TOOL_DEPENDENCIES + "/picard/2.7.0/picard.jar"
TOOL_GATK=TOOL_DEPENDENCIES + "/GATK/3.5/"
TOOL_STAR=TOOL_DEPENDENCIES + "/STAR/2.4.2a/bin/Linux_x86_64/"

TOOL_COLLECT=PROJECT_DEPENDENCIES + "/collectBpipeLogs/" // this is provided by NGSpipe2go and hence locates to the project folder
