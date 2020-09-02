// This file defines the preambles to be loaded into modules
// and the function providing the preamble to modules.
//
// How does the preambles system work:
//   * there's a default preamble defined as a global var (`default_preamble`),
//     propagated into the `module_preambles.default` key Map.
//   * modules can have custom preambles defined in the module_preambles Map. The
//     custom preable _overrides_ the default (unless it is defined as `default_preamble + custom_preamble`).
// See example in the next paragraph for more details
default_preamble="""
    export TMP="$TMP";
    [[ -n \$TMP && ! -d \$TMP ]] && mkdir -p "\$TMP";
    [[ -n \$SLURM_JOB_ID ]] && export TMP="/jobdir/\$SLURM_JOB_ID";
"""

// This map defines the default and module-specific preambles.
// The structure looks like (adding a custom preamble for bowtie as example):
//    module_preambles=[
//        default: default_preamble,
//        "bowtie": default_preamble + " && echo Running bowtie version: \$(bowtie --version | grep version)"
//    ]
// Note that bowtie's custom preamble reuses the default (otherwise, the default is overriden)
module_preambles=[
    default: default_preamble
]

//
// get_preamble
// Given a module, return its preamble string if exists, the default otherwise.
//
String get_preamble (String module) {
    return (module_preambles.containsKey(module) ? module_preambles[module] : module_preambles.default)
}

