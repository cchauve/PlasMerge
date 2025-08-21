#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
BIN_DIR=${HOME_DIR}/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv

NB_SAMPLES=$1
ALPHA=$2
MIN_LEN=100
MAX_CALLS=1000000

DATE=`date --rfc-3339=date`
TRACE_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.${DATE}.txt
RERUN_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.${DATE}.csv

BINNING="ground_truth,gplascc,mobrecon,plasbinflow"
CLASSIFICATION="plasclass,plasgraph2,rfplasmid,mlplasmids"
MERGED="unmerged,merged"
source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate
python3.9 run_utils.py check_plaseval \
	  -d ${DATA_FILE} \
	  -n ${NB_SAMPLES} \
	  -b ${BINNING} \
	  -c ${CLASSIFICATION} \
	  -bm ${MERGED} \
	  -a ${ALPHA} \
	  -ml ${MIN_LEN} \
	  -mc ${MAX_CALLS} \
	  -od ${OUTPUT_DIR} \
	  -of ${RERUN_FILE} \
	  > ${TRACE_FILE}
