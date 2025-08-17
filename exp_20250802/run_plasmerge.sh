#!/bin/bash

# Directory containing the PlasMerge python scripts
BIN_DIR=$1
# Sample ID and assembler, recorded in lowercase
IN_SAMPLE_ID=$2
IN_ASSEMBLER=$(echo ${3} | awk '{print tolower($0)}')
# Binning method, in PlasMerge syntax (gp, gt, mob, pbf)
IN_METHOD=$4
# Gzipped GFA file for assembled sample
IN_FILE_GFA=$5
# Plasmidness scores file i PlasMerge format
IN_FILE_SCORES=$6
# Input (unmerged) bins in PlasMerge format
IN_FILE_BINS=$7
# GC content bins file in PlasMerge format
IN_FILE_GC=$8
# Directory where PlasMerge results files are written
OUT_DIR=$9
# Paths of PlasMerge results files: merging scores and merged bins
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
