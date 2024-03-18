PIPELINE="test"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"    // may need adjustment for some projects

// essential vars
PROJECT="/fsimb/imbc_home/ssayolsp/tmp/test"
LOGS=PROJECT + "/logs"
OUT=PROJECT + "/out"
TMP=PROJECT + "/tmp"

// load external things
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"
load PIPELINE_ROOT + "/config/validate_module_params.groovy"

//MAIN PIPELINE TASK
test = { 
  output.dir = OUT
  def branch_outdir = new File(output.dir).getName()

  def PREAMBLE = get_preamble(module:"test", branch:branch, branch_outdir:branch_outdir)

  transform("*.in") to (".out") {
      exec """
          >&2 echo "before preamble, logs go to the screen";
          ${PREAMBLE} &&
          >&2 echo "after preamble, logs go to the corresponding file";
          cat $input > $output;
      """
  }
}

Bpipe.run { "%.in" * [ test ] }

