#!/bin/bash

# Directory where to write figures
FIG_DIR=$1
# CSV file with PlasEval results
DATA_FILE=$2

FIG_TYPE="difference"

source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate

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
PRECISION="Precision"
RECALL="Recall"
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

STAT=${NB_BINS}","${NB_CTGS}","${LEN_CTGS}
PREFIX=${FIG_TYPE}"_"${STAT}
logger -s ${FIG_TYPE} ${STAT}
OUT_FILE=${FIG_DIR}/${PREFIX}"_aggregated.png"
python PlasEval_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT}

for BINNING in ${GT} ${GP} ${MOB} ${PBF};
do
    for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
    do
	logger -s "   " ${BINNING} ${CLASSIFICATION}
	OUT_FILE=${FIG_DIR}/${PREFIX}"_"${BINNING}"_"${CLASSIFICATION}".png"
	python PlasEval_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -b ${BINNING} -c ${CLASSIFICATION}
    done
done

STAT=${PRECISION}","${RECALL}","${F1}
for VERSION in "u" "w";
do
    PREFIX=${FIG_TYPE}"_"${VERSION}"."${STAT}
    logger -s ${FIG_TYPE} ${VERSION}"."${STAT}
    OUT_FILE=${FIG_DIR}/${PREFIX}"_aggregated.png"
    python PlasEval_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -v ${VERSION}
    
    for BINNING in ${GT} ${GP} ${MOB} ${PBF};
    do
	for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
	do
	    logger -s "   " ${BINNING} ${CLASSIFICATION}
	    OUT_FILE=${FIG_DIR}/${PREFIX}"_"${BINNING}"_"${CLASSIFICATION}".png"
	    python PlasEval_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -b ${BINNING} -c ${CLASSIFICATION} -v ${VERSION}
	done
    done
done

STAT=${DISSIM}","${CUTS}","${JOINS}","${CUTS_JOINS}","${EXTRA}","${MISSING}","${CTGS}
for VERSION in "u" "n";
do
    PREFIX=${FIG_TYPE}"_"${VERSION}"."${STAT}
    logger -s ${FIG_TYPE} ${VERSION}"."${STAT}
    OUT_FILE=${FIG_DIR}/${PREFIX}"_aggregated.png"
    python PlasEval_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -v ${VERSION}
    
    for BINNING in ${GT} ${GP} ${MOB} ${PBF};
    do
	for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
	do
	    logger -s "   " ${BINNING} ${CLASSIFICATION}
	    OUT_FILE=${FIG_DIR}/${PREFIX}"_"${BINNING}"_"${CLASSIFICATION}".png"
	    python PlasEval_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -b ${BINNING} -c ${CLASSIFICATION} -v ${VERSION}
	done
    done
done
