PROJECT=$1
LOGS=$2

for f in ${PROJECT}/.bpipe/outputs/*; 
do 
 echo $f
 F=$(basename $f)
 #echo "FILEname " $F
 JOB=$(echo $F | cut -d. -f1)
 #echo "JOB " ${JOB}
 ID=$(grep -E "^commandId" $f | cut -d= -f2)
 #echo "ID " $ID
 FILE=$(grep -E "^outputFile" $f|cut -d= -f2) 
 FILE=$(basename ${FILE})                       
 echo "JOB: ${JOB}, ID: ${ID}, FILE: ${FILE}, DIR: ${PROJECT}"
 if [ ! -d "${LOGS}/${JOB}" ]; 
 then
  echo "mkdir ${LOGS}/${JOB}"
  mkdir -p ${LOGS}/${JOB}
 fi

 if [ -e ${PROJECT}/.bpipe/commandtmp/${ID}/cmd.err ]; 
 then
  echo "${PROJECT}/.bpipe/commandtmp/${ID}/cmd.err --> ${LOGS}/${JOB}/${FILE}.log"
  cp -v ${PROJECT}/.bpipe/commandtmp/${ID}/cmd.err ${LOGS}/${JOB}/${FILE}.log
  cp -v ${PROJECT}/.bpipe/commandtmp/${ID}/cmd.out ${LOGS}/${JOB}/${FILE}.out
 fi
done 
