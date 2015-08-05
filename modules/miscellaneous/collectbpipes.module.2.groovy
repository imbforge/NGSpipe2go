collectBpipeLogs = {
	doc title: "collectBpipeLogs",
		desc:  "Copy the task logs from the .bpipe dir into the logs dir",
		constraints: "",
		author: "Sergi Sayols"
	
	exec """
		${TOOL_COLLECT}/collectBpipeLogs.sh $PROJECT $LOGS
	""","collectBpipeLogs"
}

