#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
BIN_DIR=${HOME_DIR}/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv

NB_SAMPLES=$1

DATE=`date --rfc-3339=date`
TRACE_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.${DATE}.txt
RERUN_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.${DATE}.csv

BINNING="ground_truth,gplascc,mobrecon,plasbinflow"
CLASSIFICATION="plasclass,plasgraph2,rfplasmid,mlplasmids"
source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate
python3.9 run_utils.py check_plasmerge \
       -d ${DATA_FILE} \
       -n ${NB_SAMPLES} \
       -b ${BINNING} \
       -c ${CLASSIFICATION} \
       -od ${OUTPUT_DIR} \
       -of ${RERUN_FILE} \
       > ${TRACE_FILE}
