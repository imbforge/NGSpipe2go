load PIPELINE_ROOT + "config/tools.groovy"

// Tools versions and run environments
// Names should match tools.keys() in ${PIPELINE_ROOT}/config/tools.groovy
// Notes:
//   * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.
//   * Help others finding the tools here:
//       * keep the keys *sorted*
//       * use always *lowercase* if possible
tool_versions = [
    R         = [ runenv = "lmod", version= "3.6.0"        ],
    bamqc     = [ runenv = "lmod", version= "0.1.25_devel" ],
    bedtools  = [ runenv = "lmod", version= "2.27"         ],
    bowtie2   = [ runenv = "lmod", version= "2.3.4"        ],
    bowtie    = [ runenv = "lmod", version= "1.2.2"        ],
    deeptools = [ runenv = "lmod", version= "3.1"          ],
    fastqc    = [ runenv = "lmod", version= "0.11.8"       ],
    java      = [ runenv = "lmod", version= "1.8"          ],
    ketnutils = [ runenv = "lmod", version= "v365"         ],
    macs2     = [ runenv = "lmod", version= "2.1.2"        ],
    picard    = [ runenv = "lmod", version= "2.17.6"       ],
    rseqc     = [ runenv = "lmod", version= "3.0.0"        ],
    samtools  = [ runenv = "lmod", version= "1.9"          ]
]
