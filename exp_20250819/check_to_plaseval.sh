#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
BIN_DIR=${HOME_DIR}/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv

NB_SAMPLES=$1
RUN=$2
LOG_DIR=$3

DATE=`date --rfc-3339=date`
TRACE_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.convert.${RUN}.${DATE}.txt
rm -f ${TRACE_FILE}

BINNING="ground_truth,gplascc,mobrecon,plasbinflow"
CLASSIFICATION="plasclass,plasgraph2,rfplasmid,mlplasmids"
MERGED="unmerged,merged"

for B in "ground_truth" "gplascc" "mobrecon" "plasbinflow";
do
    for C in "plasclass" "plasgraph2" "rfplasmid" "mlplasmids";
    do	
	for M in "unmerged" "merged";
	do
	    logger -s "Checking conversion to PlasEval format ${B} ${C} ${M}"
	    grep ERROR ${LOG_DIR}/convert_${B}_${C}_${M}*.out >> ${TRACE_FILE}
	done
    done
done
