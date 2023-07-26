// This file defines data structures with all tools + versions + running_environments
// available to all pipelines
//
// If you install a new tool to be used in a pipeline:
//   * add the corresponding entry in tools_envs for that tool + version + runenv
//   * add default version + runenv in tools_defaults
//   * if needed, override the defaults in ${PIPELINE_ROOT}/pipelines/<pipeline>/tools.groovy
//
// Tips:
//   * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.
//   * Help others finding the tools here:
//       * keep the keys *sorted*
//       * use always *lowercase* if possible
//   * One can use only major versions in the 2nd level, and the version_str will
//     take care to report the right version. This is specially useful for programs
//     that ship often new versions with irrelevant changes.
//     For example, java 1.8: we have installed v1.8.0_101 in the lmod, and v1.8.0_181
//     in singularity. The pipeline will run version_str to report the right version,
//     but to keep this Map simple we have both under java->"1.8".
//
def conda_tools       = "/fsimb/common/conda_tools"
def conda_call        = "module try-load conda &&"  
def singularity_tools = "/fsimb/common/singularity_tools"

// default runenvs and versions for each tools.
// Names should match those of tools_envs

tools_defaults = [
    R              : [ runenv: "lmod", version: "R/Bioconductor_3.14"],
    bamqc          : [ runenv: "lmod", version: "0.1.25"       ],
    bamutil        : [ runenv: "lmod", version: "1.0.14"       ],
    bedtools       : [ runenv: "lmod", version: "2.29.2"       ],
    bowtie         : [ runenv: "lmod", version: "1.3.1"        ],
    bowtie2        : [ runenv: "lmod", version: "2.4.5"        ],
    cellranger     : [ runenv: "lmod", version: "6.0.0"        ],
    cellrangerarc  : [ runenv: "lmod", version: "2.0.0"        ],
    cellrangeratac : [ runenv: "lmod", version: "2.0.0"        ],
    bwa            : [ runenv: "lmod", version: "0.7.17"       ],
    cutadapt       : [ runenv: "lmod", version: "4.0"          ],
    deeptools      : [ runenv: "lmod", version: "3.5.1_singularity"],
    fastqc         : [ runenv: "lmod", version: "0.11.9"       ],
    fastqscreen    : [ runenv: "lmod", version: "0.15"         ],
    fastx          : [ runenv: "lmod", version: "0.0.14"       ],
    gatk           : [ runenv: "lmod", version: "4.2.6.1"      ],
    htseq          : [ runenv: "lmod", version: "2.0.2"        ],
    java           : [ runenv: "lmod", version: "1.8"          ],
    kentutils      : [ runenv: "lmod", version: "v385"         ],
    macs2          : [ runenv: "lmod", version: "2.1.2"        ],
    mirdeep2       : [ runenv: "lmod", version: "2.0.1.3"      ],
    multiqc        : [ runenv: "lmod", version: "1.9"          ],
    pear           : [ runenv: "lmod", version: "0.9.11"       ],
    picard         : [ runenv: "lmod", version: "2.20"         ],
    pingpongpro    : [ runenv: "lmod", version: "1.0"          ],
    qualimap       : [ runenv: "lmod", version: "2.2.1"        ],
    repenrich      : [ runenv: "lmod", version: "1.2"          ],
    rmats          : [ runenv: "lmod", version: "4.1.2"        ],
    rseqc          : [ runenv: "lmod", version: "4.0.0"        ],
    samtools       : [ runenv: "lmod", version: "1.10"         ],
    seqtk          : [ runenv: "lmod", version: "1.3"          ], 
    starfusion     : [ runenv: "lmod", version: "0.8.0"        ],
    star           : [ runenv: "lmod", version: "2.7"          ],
    stringtie      : [ runenv: "lmod", version: "2.2.1"        ],
    subread        : [ runenv: "lmod", version: "2.0.0"        ],
    trimgalore     : [ runenv: "lmod", version: "0.5.0"        ],
    umitools       : [ runenv: "lmod", version: "1.1.2"        ]
]





// This map defines how to prepare the environment in order to have in PATH all
// tools + versions + running_environments available. Structure:
//   * 1st level: tool name
//   * 2nd level: version --> version *string*
//   * 3rd level: runenv  --> one of [lmod|singularity|conda]

tools_envs = [
    R: [
        "3.6.0": [
            singularity: "alias Rscript=\"singularity run --app Rscript ${singularity_tools}/R/3.6.0r0/R.simg\""
        ],
        "4.0.3": [
            lmod: "module load R/4.0.3"
        ],
        "4.2.0": [
            lmod: "module load R/4.2.0"
        ],
        "R/Bioconductor_3.14" : [
            lmod: "module load R/Bioconductor_3.14_singularity"
        ]
    ],
    bamqc: [
        "0.1.25": [
            lmod: "module load BamQC/0.1.25_devel"
        ]
    ],
    bamutil: [
        "1.0.13": [
            lmod: "module load bamUtil/1.0.13"
        ],
        "1.0.14": [
            lmod: "module load bamUtil/1.0.14"
        ]
    ],
    bedtools: [
        "2.27": [
            lmod: "module load bedtools/2.27.1",
            conda: "source activate ${conda_tools}/bedtools/2.27.1",
            singularity: "alias bedtools=\"singularity run --app bedtools ${singularity_tools}/bedtools/2.27.1r0/bedtools.simg\""
        ],
        "2.29": [
            lmod: "module load bedtools/2.29.2"
        ]
    ],
    bowtie: [
        "1.2.2": [
            lmod: "module load bowtie/1.2.2",
            conda: "source activate ${conda_tools}/bowtie/1.2.2"
        ],
        "1.3.1": [
            lmod: "module load bowtie/1.3.1"
        ]
    ],
    bowtie2: [
        "2.3.4": [
            lmod: "module load bowtie2/2.3.4.3",
            conda: "source activate ${conda_tools}/bowtie2/2.3.4",
            singularity: "alias bowtie2=\"singularity run --app bowtie2 ${singularity_tools}/bowtie2/2.3.4.3r0/bowtie2.simg\""
        ],
        "2.4.5": [
            lmod: "module load bowtie2/2.4.5"
        ]
    ],
    bwa: [
        "0.7.17": [
            lmod: "module load bwa/0.7.17"
        ]
    ],
    cellranger: [
        "6.0.0": [
            lmod: "module load cellranger/6.0.0"
        ],
        "7.0.1": [
            lmod: "module load cellranger/7.0.1"
        ]
    ],
    cellrangerarc: [
        "2.0.0": [
            lmod: "module load cellrangerARC/2.0.0"
        ]
    ],
    cellrangeratac: [
        "2.0.0": [
            lmod: "module load cellrangerATAC/2.0.0"
        ]
    ],
    cutadapt: [
        "1.18": [
            conda: "source activate ${conda_tools}/cutadapt/1.18",
            singularity :"alias cutadapt=\"singularity run --app cutadapt ${singularity_tools}/cutadapt/1.18r0/cutadapt.simg\""
        ],
        "4.0": [
            lmod: "module load cutadapt/4.0"
        ],
        "4.4": [
            lmod: "module load cutadapt/4.4_singularity"
        ]
    ],
    deeptools: [
        "3.5.1": [
            lmod: "module load deepTools/3.5.1_singularity"
        ]
    ],
    fastqc: [
        "0.11.8": [
            conda: "source activate ${conda_tools}/fastqc/0.11.8",
            singularity: "alias fastqc=\"singularity run --app fastqc ${singularity_tools}/fastqc/0.11.8r0/fastqc.simg\""
        ],
        "0.11.9": [
            lmod: "module load fastqc/0.11.9"
        ]
    ],
    fastqscreen: [
        "0.13": [
            lmod: "module load fastq_screen/0.13.0",
            conda: "source activate ${conda_tools}/fastq_screen/0.13"
        ],
        "0.14": [
            lmod: "module load fastq_screen/0.14.1"
        ],
        "0.15": [
            lmod: "module load fastq_screen/0.15.2"
        ]
    ],
    fastx: [
        "0.0.14": [
            lmod: "module load fastx_toolkit/0.0.14",
            conda: "source activate ${conda_tools}/fastx_toolkit/0.0.14"
        ]
    ],
    gatk: [
        "4.2.6.1": [
            lmod: "module load GATK/4.2.6.1_singularity"
        ]
    ],
    htseq: [
        "0.6.1": [
            lmod: "module load htseq/0.6.1"
        ]
    ],
    java: [
        "1.8": [
            lmod: "module load jdk/1.8.0_332",
            singularity: "alias java=\"singularity run --app java ${singularity_tools}/openjdk/8u181r0/openjdk8.simg\""
        ]
    ],
    kentutils: [
        "v385": [
            lmod: "module load kentUtils/v385"
        ]
    ],
    macs2: [
        "2.1.2": [
            lmod: "module load macs2/2.1.2",
            conda: "source activate ${conda_tools}/macs2/2.1.2",
            singularity: "alias macs2=\"singularity run --app macs2 ${singularity_tools}/macs2/2.1.2.1r0/macs.simg\""
        ]
    ],
    mirdeep2: [
        "2.0.1.3": [
            lmod: "module load mirdeep2/2.0.1.3",
            conda: "source activate ${conda_tools}/mirdeep2/2.0.1.3"
        ]
    ],
    multiqc: [
        "1.7": [
            lmod: "module load MultiQC/1.7",
            conda: "${conda_call} source activate ${conda_tools}/MultiQC/1.7",
            singularity: "alias multiqc=\"singularity run --app multiqc ${singularity_tools}/MultiQC/1.7/multiqc.simg\""
        ],
        "1.9": [
            lmod: "module load MultiQC/1.9"
        ]
    ],
    pear: [
        "0.9.11": [
            lmod: "module load pear/0.9.11"
        ]
    ],
    picard: [
        "2.7": [
            lmod: "module load picard/2.7.0"
        ],
        "2.18": [
            conda: "source activate ${conda_tools}/picard/2.18.26",
            singularity: "alias picard=\"singularity run --app picard ${singularity_tools}/picard/2.18.17r0/picard.simg\""
        ],
        "2.20": [
            lmod: "module load picard/2.20.3"
        ]
    ],
    pingpongpro: [
        "1.0": [
            lmod: "module load pingpongpro/1.0"
        ]
    ],
    qualimap: [
        "2.2.1": [
            lmod: "module load qualimap/2.2.1"
        ]
    ],
    repenrich: [
        "1.2": [
            conda: "source activate ${conda_tools}/repenrich/1.2"
        ]
    ],
    rmats: [
        "4.1.2": [
            lmod: "module load rmats/4.1.2"
        ]
    ],
    rseqc: [
        "3.0.0": [
            lmod: "module load RSeQC/3.0.0",
            conda: "source activate ${conda_tools}/rseqc/3.0.0"
        ],
        "4.0.0": [
            lmod: "module load RSeQC/4.0.0"        
        ]
    ],
    samtools: [
        "1.9": [
            lmod: "module load samtools/1.9",
            conda: "source activate ${conda_tools}/samtools/1.9",
            singularity: "alias samtools=\"singularity run --app samtools ${singularity_tools}/samtools/1.9r1/samtools.simg\""
        ],
        "1.10": [
            lmod: "module load samtools/1.10"
        ]
    ],
    seqtk: [
        "1.3": [
            lmod: "module load seqtk/1.3",
            conda: "source activate ${conda_tools}/seqtk/1.3"
        ]
    ],
    starfusion: [
        "0.8.0": [
            lmod: "module load STAR-Fusion/0.8.0"
        ]
    ],
    star: [
        "2.7.3": [
            lmod: "module load star/2.7.3a",
            conda: "source activate ${conda_tools}/star/2.7.0d",
            singularity: "alias STAR=\"singularity run --app STAR ${singularity_tools}/star/2.7.0fr0/STAR.simg\""
        ],
        "2.7.10": [
            lmod: "module load star/2.7.10a"
        ]
    ],
    stringtie: [
        "2.2.1": [
            lmod: "module load stringtie/2.2.1"
        ]
    ],
    subread: [
        "1.6": [
            conda: "source activate ${conda_tools}/subread/1.6.3",
            singularity: """
                alias exactSNP=\"singularity run --app exactSNP ${singularity_tools}/subread/1.6.3r0/subread.simg\";
                alias featureCounts=\"singularity run --app featureCounts ${singularity_tools}/subread/1.6.3r0/subread.simg\";
                alias subindel=\"singularity run --app subindel ${singularity_tools}/subread/1.6.3r0/subread.simg\";
                alias subjunc=\"singularity run --app subjunc ${singularity_tools}/subread/1.6.3r0/subread.simg\";
                alias sublong=\"singularity run --app sublong ${singularity_tools}/subread/1.6.3r0/subread.simg\";
                alias subread-align=\"singularity run --app subread-align ${singularity_tools}/subread/1.6.3r0/subread.simg\";
                alias subread-buildindex=\"singularity run --app subread-buildindex${singularity_tools}/subread/1.6.3r0/subread.simg\"
            """
        ],
        "2.0": [
            lmod: "module load subread/2.0.0"
        ]
    ],
    trimgalore: [
        "0.5.0": [
            singularity: "alias trim_galore=\"singularity run --app trim_galore ${singularity_tools}/trimgalore/0.5.0r0/trimgalore.simg\""
        ]
    ],
    umitools: [
        "1.1.2": [
            lmod: "module load umitools/1.1.2_singularity"
        ]
    ]
]

//
// prepare_tool_env
// 
// Given a tool, version and run env, return its tools_envs string.
// If any of the keys doesn't exist, returns the POSIX shell null command.
//
String prepare_tool_env (String tool, String version, String runenv) {

    if(tools_envs.containsKey(tool) &&
       tools_envs[tool].containsKey(version) &&
       tools_envs[tool][version].containsKey(runenv)) {
        return tools_envs[tool][version][runenv]
    } else {
        return ":"   // return the POSIX shell null command
    }
}


