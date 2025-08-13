#!/bin/bash

# source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate
# module load StdEnv/2020 gurobi/9.1.2

BIN_DIR=$1
IN_SAMPLE_ID=$2
IN_ASSEMBLER=$(echo ${3} | awk '{print tolower($0)}')
IN_METHOD=$4
IN_FILE_GFA=$5
IN_FILE_SCORES=$6
IN_FILE_BINS=$7
IN_FILE_GC=$8
OUT_DIR=$9
OUT_FILE_SCORES=${10}
OUT_FILE_BINS=${11}

EXP_ID=${IN_SAMPLE_ID}.${IN_ASSEMBLER}

mkdir -p ${OUT_DIR}

EXIT=0

for FILE in ${IN_FILE_GFA} ${IN_FILE_SCORES} ${IN_FILE_BINS} ${IN_FILE_GC}
do
    if ! [ -f ${FILE} ];
    then
	logger -s "ERROR.MISSING" ${EXP_ID} ${FILE}
	EXIT=1
    else
	FILE_SIZE=$(stat -c%s "${FILE}")
	if (( ${FILE_SIZE} == 0 ));
	then
	    logger -s "ERROR.EMPTY" ${EXP_ID} ${FILE}
	    EXIT=1
	fi
    fi
done

if (( ${EXIT} == 0 ));
then
    python3.9 ${BIN_DIR}/plasmerge.py \
 	      -s ${IN_SAMPLE_ID} \
 	      -a ${IN_FILE_GFA} \
	      -t ${IN_ASSEMBLER} \
 	      -p ${IN_FILE_SCORES} \
 	      -b ${IN_FILE_BINS} \
 	      -r ${IN_METHOD} \
 	      -d ${OUT_DIR} \
 	      -l \
 	      -os ${OUT_FILE_SCORES} \
	      -om ${OUT_FILE_BINS} \
 	      -g ${IN_FILE_GC}
    
    ERROR=0

    for FILE in ${OUT_FILE_SCORES} ${OUT_FILE_BINS}
    do
	if ! [ -f ${FILE} ];
	then
	    logger -s "ERROR.MISSING" ${EXP_ID} ${FILE}
	    ERROR=1
	else
	    FILE_SIZE=$(stat -c%s "${FILE}")
	    if (( ${FILE_SIZE} == 0 ));
	    then
		logger -s "ERROR.EMPTY" ${EXP_ID} ${FILE}
		ERROR=1
	    fi
	fi
    done
    
    if (( ${ERROR} == 0 ));
    then
	logger -s "LOG.OUTPUT" ${EXP_ID} ${OUT_FILE_SCORES} ${OUT_FILE_BINS}
    fi
else
    logger -s "ERROR.INPUT" ${EXP_ID} "missing or empty data file"
fi
