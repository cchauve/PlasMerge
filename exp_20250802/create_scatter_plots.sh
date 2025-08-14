#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge/
BIN_DIR=/${HOME_DIR}/src/
EXP_DIR=${HOME_DIR}/exp_20250802/
#INPUT_DIR=${EXP_DIR}/eval/
INPUT_DIR=${EXP_DIR}/eval_0/
OUTPUT_DIR=${EXP_DIR}/analysis/
#FIG_DIR=${OUTPUT_DIR}/figures/
FIG_DIR=${OUTPUT_DIR}/figures_0/
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
CUTS_JOINS="Cuts,Joins"
EXTRA="Extra_ctgs"
MISSING="Missing_ctgs"
CTGS="Extra_ctgs,Missing_ctgs"
NB_BINS="nb_bins"
NB_CTGS="nb_ctgs"
LEN_CTGS="len_ctgs"

source ${HOME_DIR}/../env_plaseval/bin/activate

#for STAT in ${PREC} ${REC} ${F1} ${DISSIM} ${CUTS} ${JOINS} ${CUTS_JOINS} ${EXTRA} ${MISSING} ${CTGS} ${NB_BINS} ${NB_CTGS} ${LEN_CTGS}
for STAT in ${DISSIM} ${CUTS} ${JOINS} ${CUTS_JOINS} ${EXTRA} ${MISSING} ${CTGS} 
do
    logger -s ${STAT}
    python analysis_utils.py \
	   scatter_aggregated \
	   ${DATA_FILE} \
	   ${FIG_DIR} \
	   ${STAT}

    for BINNING in ${GT} ${GP} ${MOB} ${PBF};
    do
	for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
	do
	    logger -s "   " ${BINNING} ${CLASSIFICATION}
	    python analysis_utils.py \
		   scatter_combination \
		   ${DATA_FILE} \
		   ${FIG_DIR} \
		   ${STAT} \
		   ${BINNING} \
		   ${CLASSIFICATION}
	done
    done
done
