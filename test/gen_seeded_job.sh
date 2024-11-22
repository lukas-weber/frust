#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: $0 JOBFILE SEEDED_JOBFILE
	exit 1
fi
sed $1 -e 's/tm\.sweeps.*$/tm.sweeps = 50\ntm.seed = 18088206/g' \
       -e 's/tm\.thermalization.*$/tm.thermalization = 50/g' \
       -e 's/tm\.binsize.*$/tm.binsize = 1/g' > $2
chmod +x $2
