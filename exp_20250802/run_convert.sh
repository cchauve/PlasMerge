#!/bin/bash

# source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate
# module load StdEnv/2020 gurobi/9.1.2

BIN_DIR=$1
IN_FILE=$2
OUT_FILE=$3
IN_METHOD=$4

if [ -f ${IN_FILE} ]; then
    logger -s "### Input file " ${IN_FILE}
    python3.9 ${BIN_DIR}/convert_utils.py ${IN_FILE} ${OUT_FILE} ${IN_METHOD}
    if ! [ -f ${OUT_FILE} ]; then
	logger -s "ERROR no output file " ${OUT_FILE}
    else
	logger -s "SUCCESS output file " ${OUT_FILE}
    fi
else
    logger -s "ERROR missing input file " ${IN_FILE}
fi
