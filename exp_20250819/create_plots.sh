#!/bin/bash

# Stats to plot: bins, eval, comp
MODE=$1
# Figure type: "scatter" or "difference"
FIG_TYPE=$2
# Directory where to write figures
FIG_DIR=$3
# CSV file with PlasEval results
DATA_FILE=$4
# Binning methods
BINNING_METHODS=$5
# Classification methods
CLASSIFICATION_METHODS=$6
# Alpha value, used only for comp, "NA" otherwise
ALPHA=$7
# Type of stats: "u", "w" or "u w" for eval, "u", "n" or "u n" for comp, "NA" for bins
NORM=$8

source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate

# eval stats
PRECISION="Precision"
RECALL="Recall"
F1="F1"
# comp stats
DISSIM="Dissimilarity"
CUTS="Cuts"
JOINS="Joins"
CUTS_JOINS="Cuts,Joins"
EXTRA="Extra_ctgs"
MISSING="Missing_ctgs"
CTGS="Extra_ctgs,Missing_ctgs"
# Bins stats
NB_BINS="nb_bins"
NB_CTGS="nb_ctgs"
LEN_CTGS="len_ctgs"

if [ "${MODE}" = "bins" ]; then
    for STAT in ${NB_BINS} ${NB_CTGS} ${LEN_CTGS};
    do
	PREFIX=${FIG_TYPE}"_"`echo ${STAT} | sed 's/,/_/g'`
	OUT_FILE=${FIG_DIR}/${PREFIX}"_aggregated.png"
	logger -s ${FIG_TYPE} ${STAT} ${OUT_FILE} 
	python analysis_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT}

	for BINNING in ${BINNING_METHODS};
	do
	    for CLASSIFICATION in ${CLASSIFICATION_METHODS};
	    do
		OUT_FILE=${FIG_DIR}/${PREFIX}"_"${BINNING}"_"${CLASSIFICATION}".png"
		logger -s "   " ${BINNING} ${CLASSIFICATION} ${OUT_FILE} 
		python analysis_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -b ${BINNING} -c ${CLASSIFICATION}
	    done
	done
    done
fi

if [ "${MODE}" = "eval" ]; then
    for STAT in ${PRECISION} ${RECALL} ${F1};
    do
	for VERSION in ${NORM};
	do
	    PREFIX=${FIG_TYPE}"_"${VERSION}"."`echo ${STAT} | sed 's/,/_/g'`
	    logger -s ${FIG_TYPE} ${VERSION}"."${STAT}
	    OUT_FILE=${FIG_DIR}/${PREFIX}"_aggregated.png"
	    python analysis_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -v ${VERSION}
	    
	    for BINNING in ${BINNING_METHODS};
	    do
		for CLASSIFICATION in ${CLASSIFICATION_METHODS};
		do
		    OUT_FILE=${FIG_DIR}/${PREFIX}"_"${BINNING}"_"${CLASSIFICATION}".png"
		    logger -s "   " ${BINNING} ${CLASSIFICATION} ${OUT_FILE}
		    python analysis_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -b ${BINNING} -c ${CLASSIFICATION} -v ${VERSION}
		done
	    done
	done
    done
fi

if [ "${MODE}" = "comp" ]; then
    for STAT in ${DISSIM} ${CUTS} ${JOINS} ${CUTS_JOINS} ${EXTRA} ${MISSING} ${CTGS};
    do
	for VERSION in ${NORM};
	do
	    PREFIX=${FIG_TYPE}"_"${VERSION}"."`echo ${STAT} | sed 's/,/_/g'`"."${ALPHA}
	    OUT_FILE=${FIG_DIR}/${PREFIX}"_aggregated.png"
	    logger -s ${FIG_TYPE} ${VERSION}"."${STAT} ${OUT_FILE}
	    python analysis_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -v ${VERSION}
	    
	    for BINNING in ${BINNING_METHODS};
	    do
		for CLASSIFICATION in ${CLASSIFICATION_METHODS};
		do
		    OUT_FILE=${FIG_DIR}/${PREFIX}"_"${BINNING}"_"${CLASSIFICATION}".png"
		    logger -s "   " ${BINNING} ${CLASSIFICATION} ${OUT_FILE}
		    python analysis_utils.py ${FIG_TYPE} ${DATA_FILE} ${OUT_FILE} ${STAT} -b ${BINNING} -c ${CLASSIFICATION} -v ${VERSION}
		done
	    done
	done
    done
fi
