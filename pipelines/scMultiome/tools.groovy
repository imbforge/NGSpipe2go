// Tools custom versions and run environments
// Overrides the defaults defined in ${PIPELINE_ROOT}/config/tools.groovy
// Names should match tools_defaults.keys() in ${PIPELINE_ROOT}/config/tools.groovy
// 
// The structure of this map is:
//   tools_custom = [
//       R       : [ runenv: "lmod", version: "3.6.0" ],
//    <...>
//       samtools: [ runenv: "lmod", version: "1.9"   ]
//   ]
//
// Tips:
//   * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.
load PIPELINE_ROOT + "/config/tools.groovy" // tools_defaults are specified here
tools_custom = [ ] 

tools = new LinkedHashMap(tools_defaults)   // create new tools map based on defaults
tools.putAll(tools_custom)                  // override with users custom versions/runenvs
