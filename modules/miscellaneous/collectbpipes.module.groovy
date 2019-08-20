collectBpipeLogs = {
	doc title: "collectBpipeLogs",
		desc:  "Copy the task logs from the .bpipe dir into the logs dir",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"
	
    def PREAMBLE = get_preamble("collectBpipeLogs")

	exec """
        ${PREAMBLE} &&

		for f in $PROJECT/.bpipe/outputs/*; do
			F=\$(basename \$f)                              &&
			JOB=\$(echo \$F |cut -d. -f1)                   &&
			ID=\$(grep -E "^commandId" \$f|cut -d= -f2)     &&
			FILE=\$(grep -E "^outputFile" \$f|cut -d= -f2)  &&
			FILE=\$(basename \$FILE)                        &&
			echo "JOB: \${JOB}, ID: \${ID}, FILE: \${FILE}" &&
			if [ ! -d "$LOGS/\${JOB}" ]; then
				echo "mkdir $LOGS/\${JOB}"                  &&
				mkdir -p $LOGS/\${JOB}                      ;
			fi                                              ;
			if [ -e $PROJECT/.bpipe/commandtmp/\${ID}/cmd.err ]; then
				echo "$PROJECT/.bpipe/commandtmp/\${ID}/cmd.err --> $LOGS/\${JOB}/\${FILE}.log" &&
				cp $PROJECT/.bpipe/commandtmp/\${ID}/cmd.err $LOGS/\${JOB}/\${FILE}.log ;
			fi ;
		done
	""","collectBpipeLogs"
}

