// This file defines the preambles to be loaded into modules
// and the function providing the preamble to modules.
//
// How does the preambles system work:
//   * there's a default preamble defined as a global var (`default_preamble`),
//     propagated into the `module_preambles.default` key Map.
//   * modules can have custom preambles defined in the module_preambles Map. The
//     custom preable _overrides_ the default (unless it is defined as `default_preamble + custom_preamble`).
// See example in the next paragraph for more details
//
// Preambles can have placeholder variables, named within double underscore(eg. __branch__)
// placeholder variables get replaced at runtime with the actual value. The
// replacement is done in `get_preamble()`
default_preamble="""
    export TMP="$TMP";
    [[ -n \$TMP && ! -d \$TMP ]] && mkdir -p "\$TMP";
    [[ -n \$SLURM_JOB_ID ]] && export TMP="/jobdir/\$SLURM_JOB_ID";

    export LOGS="${LOGS}/__stage__/";
    [[ -n \$LOGS && ! -d \$LOGS ]] && mkdir -p "\$LOGS";
    pref=\$([[ -z "__outdir__" ]] || realpath --relative-to=. "__outdir__" | sed 's/\\//_/g');
    readonly LOG_FILE=\${LOGS}/__input___\${pref}.log;
    touch \$LOG_FILE;
    exec 1>\$LOG_FILE;
    exec 2>&1;
    :
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
// Args: the function receives a Map, with at least 2 elements:
//   stage: the module that is calling the function
//   outdir: the output directory
//   input: the basename of the input file name
// Other args are optional, and will be used to replace the placeholder variables
// used within the preamble.
//
String get_preamble(Map args) {

  // fetch the preamble for the module, if this exists. Otherwise, pick the default preamble
  def x = (module_preambles.containsKey(args.module) ? module_preambles[args.module] : module_preambles.default)

  // replace placeholder variables within the preamble with the actual runtime value
  args.each { key, value ->
    x = x.replaceAll("__" + key + "__", value)
  }

  return(x)
}


