#!/bin/bash

HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
EXP_DIR=${HOME_DIR}/exp_20250819
RESULTS_DIR=${EXP_DIR}/output
ANALYSIS_DIR=${EXP_DIR}/analysis
FIG_DIR=${EXP_DIR}/figures
mkdir -p ${ANALYSIS_DIR}
mkdir -p ${FIG_DIR} 
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv

source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate

NB_SAMPLES=$1
RUN=$2
ALPHA=$3
# if ALPHA is eval, only PlasEval eval results are read

DATE=`date --rfc-3339=date`
EVAL_FILE=${ANALYSIS_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.${ALPHA}.${RUN}.${DATE}
EVAL_FILE_CSV=${EVAL_FILE}.csv
EVAL_FILE_NA=${EVAL_FILE}.NA.txt

# logger -s "Creating CSV file for run1 without PlasEval comp results"
# python analysis_utils.py csv \
#        ${DATA_FILE} \
#        ${RESULTS_DIR} \
#        ${EVAL_FILE_CSV} \
#        -n ${NB_SAMPLES} \
#        -a ${ALPHA} \
#        -v \
#        > ${EVAL_FILE_NA}

BINNING="ground_truth gplascc mobrecon plasbinflow"
CLASSIFICATION="plasclass plasgraph2 mlplasmids rfplasmid"

logger -s "Creating plots for bins statistics"
./create_plots.sh bins scatter ${FIG_DIR} ${EVAL_FILE_CSV} "${BINNING}" "${CLASSIFICATION}" "NA" "NA"
./create_plots.sh bins difference ${FIG_DIR} ${EVAL_FILE_CSV} "${BINNING}" "${CLASSIFICATION}" "NA" "NA"

logger -s "Creating plots for eval statistics"
./create_plots.sh eval scatter ${FIG_DIR} ${EVAL_FILE_CSV} "${BINNING}" "${CLASSIFICATION}" "NA" "u w"
./create_plots.sh eval difference ${FIG_DIR} ${EVAL_FILE_CSV} "${BINNING}" "${CLASSIFICATION}" "NA" "u w"
