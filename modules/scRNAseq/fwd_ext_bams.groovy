fwd_ext_bams = {
        doc title: "Forward external bam files",
        desc:  "Find bam files within an external result directory and forward to nex module",
        constraints: "Designed for result folder structure from ScaleRNA software from ScaleBio ",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

EXTERNAL_RESULTDIR = input // must be global because used by other modules

    def static_subdir = (fwd_ext_bams_vars.static_subdir  ? fwd_ext_bams_vars.static_subdir : "")
    def wildcard_subdir = (fwd_ext_bams_vars.static_subdir  ? fwd_ext_bams_vars.wildcard_subdir : "")

    def base = new File(input + "/" + fwd_ext_bams_vars.static_subdir + "/")
    def bams = []
    
    base.eachFileRecurse { file ->
        if (file.name.endsWith(".bam") && file.path.contains(fwd_ext_bams_vars.wildcard_subdir)) {
            bams << file.absolutePath
        }
    }

	println "bam files forwarded from $EXTERNAL_RESULTDIR directory: $bams"
	forward bams
}

fwd_ext_dir = { 
    println "Process results from external directory: $input"
    forward input
}
