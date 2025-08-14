#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge/
BIN_DIR=/${HOME_DIR}/src/
EXP_DIR=${HOME_DIR}/exp_20250802/
INPUT_DIR=${EXP_DIR}/eval/
OUTPUT_DIR=${EXP_DIR}/analysis/
FIG_DIR=${OUTPUT_DIR}/figures/
DATA_FILE=${OUTPUT_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv

# Methods
RFPLASMID="rfplasmid"
PLASCLASS="plasclass"
PLASGRAPH="plasgraph2"
MLPLASMIDS="mlplasmids"
GT="ground_truth"
GP="gplascc"
MOB="mobrecon"
PBF="plasbinflow"
# Statistics
PREC="Precision"
REC="Recall"
F1="F1"
DISSIM="Dissimilarity"
CUTS="Cuts"
JOINS="Joins"
EXTRA="Extra_ctgs"
MISSING="Missing_ctgs"
NB_BINS="nb_bins"
NB_CTGS="nb_ctgs"
LEN_CTGS="len_ctgs"
STATS=${PREC}","${REC}","${F1}","${DISSIM}","${CUTS}","${JOINS}","${EXTRA}","${MISSING}

source ${HOME_DIR}/../env_plaseval/bin/activate

logger -s ${STATS}
python analysis_utils.py \
       diff_aggregated \
       ${DATA_FILE} \
       ${FIG_DIR} \
       ${STATS}

for BINNING in ${GT} ${GP} ${MOB} ${PBF};
do
    for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
    do
	logger -s "   " ${BINNING} ${CLASSIFICATION}
	python analysis_utils.py \
	       diff_combination \
	       ${DATA_FILE} \
	       ${FIG_DIR} \
	       ${STATS} \
	       ${BINNING} \
	       ${CLASSIFICATION}
    done
done
