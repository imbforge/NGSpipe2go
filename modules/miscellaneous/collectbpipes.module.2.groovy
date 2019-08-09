collectBpipeLogs = {
	doc title: "collectBpipeLogs",
		desc:  "Copy the task logs from the .bpipe dir into the logs dir",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"
	
	exec """
		${PIPELINE_ROOT}/tools/collectBpipeLogs/collectBpipeLogs.sh $PROJECT $LOGS
	""","collectBpipeLogs"
}

