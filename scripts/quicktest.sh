#!/bin/bash

if [ $# -eq 0 ]; then
	echo Usage: $0 JOBFILE [ARGS...]
	exit 1
fi

ed_code=$(dirname $0)/ed.py
ed_compare=$(dirname $0)/ed_compare.py
job=$1

shift
jobname=$(basename -s .py $job)

loadl delete $job
python3 $ed_code $job -o $jobname.npz
loadl run $job
if loadl status $job --need-merge; then
    loadl merge $job
fi

python3 $ed_compare $jobname.npz $jobname.results.json $@

