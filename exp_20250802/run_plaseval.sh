#!/bin/bash

# cd /scratch/chauvec/PLASMERGE
# module load StdEnv/2020 python/3.9.6
# python -m venv  env_plaseval
# source env_plaseval/bin/activate
# pip install networkx
# pip install psutil
# pip install pandas
# pip install bidict
# pip install biopython
# source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate
# module load StdEnv/2020 gurobi/9.1.2

PLASEVAL_BIN_DIR=$1
PLASMERGE_BIN_DIR=$2
IN_SAMPLE_ID=$3
IN_ASSEMBLER=$4
IN_METHOD=$5
IN_FILE_GFA=$6
IN_FILE_BINS_PRED_UNMERGED=$7
IN_FILE_BINS_PRED_MERGED=$8
IN_FILE_BINS_GT=$9
OUT_DIR=${10}
IN_MIN_LEN=${11}
IN_ALPHA=${12}
IN_MAX_RECURSIVE_CALLS=${13}

EXP_ID=${IN_SAMPLE_ID}.${IN_ASSEMBLER}
OUT_EXP_ID=${IN_SAMPLE_ID}_${IN_ASSEMBLER}

IN_FILE_BINS_PRED_UNMERGED_REFORMATTED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.unmerged.tsv
IN_FILE_BINS_PRED_MERGED_REFORMATTED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.merged.tsv
OUT_FILE_EVAL_UNMERGED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.unmerged.eval.out
OUT_FILE_COMP_OUT_UNMERGED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.unmerged.comp.out
OUT_FILE_COMP_LOG_UNMERGED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.unmerged.comp.log
OUT_FILE_EVAL_MERGED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.merged.eval.out
OUT_FILE_COMP_OUT_MERGED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.merged.comp.out
OUT_FILE_COMP_LOG_MERGED=${OUT_DIR}/${OUT_EXP_ID}_${IN_METHOD}.merged.comp.log

logger -s "Creating output directory" ${OUT_DIR}
mkdir -p ${OUT_DIR}

EXIT=0

for FILE in ${IN_FILE_GFA} ${IN_FILE_BINS_PRED_UNMERGED} ${IN_FILE_BINS_PRED_MERGED} ${IN_FILE_BINS_GT}
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
    logger -s "Converting unmerged bins" ${IN_FILE_BINS_PRED_UNMERGED} ${IN_FILE_BINS_PRED_UNMERGED_REFORMATTED}
    python3.9 ${PLASMERGE_BIN_DIR}/convert_utils.py \
	      ${IN_FILE_BINS_PRED_UNMERGED} \
	      ${IN_FILE_BINS_PRED_UNMERGED_REFORMATTED} \
	      to_plaseval \
	      ${IN_FILE_GFA}

    logger -s "Converting merged bins" ${IN_FILE_BINS_PRED_MERGED} ${IN_FILE_BINS_PRED_MERGED_REFORMATTED}
    python3.9 ${PLASMERGE_BIN_DIR}/convert_utils.py \
	      ${IN_FILE_BINS_PRED_MERGED} \
	      ${IN_FILE_BINS_PRED_MERGED_REFORMATTED} \
	      to_plaseval \
	      ${IN_FILE_GFA}

    logger -s "Eval unmerged bins" ${OUT_FILE_EVAL_UNMERGED}
    python3.9 ${PLASEVAL_BIN_DIR}/plaseval.py \
 	      eval \
	      --pred ${IN_FILE_BINS_PRED_UNMERGED_REFORMATTED} \
	      --gt ${IN_FILE_BINS_GT} \
	      --out_file ${OUT_FILE_EVAL_UNMERGED} \
	      --min_len ${IN_MIN_LEN}
    
    logger -s "Eval merged bins" ${OUT_FILE_EVAL_MERGED}
    python3.9 ${PLASEVAL_BIN_DIR}/plaseval.py \
 	      eval \
	      --pred ${IN_FILE_BINS_PRED_MERGED_REFORMATTED} \
	      --gt ${IN_FILE_BINS_GT} \
	      --out_file ${OUT_FILE_EVAL_MERGED} \
	      --min_len ${IN_MIN_LEN}

    logger -s "Comp unmerged bins" ${OUT_FILE_COMP_OUT_UNMERGED} ${OUT_FILE_COMP_LOG_UNMERGED}
    python3.9 ${PLASEVAL_BIN_DIR}/plaseval.py \
 	      comp \
	      --l ${IN_FILE_BINS_PRED_UNMERGED_REFORMATTED} \
	      --r ${IN_FILE_BINS_GT} \
	      --out_file ${OUT_FILE_COMP_OUT_UNMERGED} \
	      --log_file ${OUT_FILE_COMP_LOG_UNMERGED} \
	      --min_len ${IN_MIN_LEN} \
	      --p ${IN_ALPHA} \
	      --max_calls ${IN_MAX_RECURSIVE_CALLS}

    logger -s "Comp merged bins" ${OUT_FILE_COMP_OUT_MERGED} ${OUT_FILE_COMP_LOG_MERGED}
    python3.9 ${PLASEVAL_BIN_DIR}/plaseval.py \
 	      comp \
	      --l ${IN_FILE_BINS_PRED_MERGED_REFORMATTED} \
	      --r ${IN_FILE_BINS_GT} \
	      --out_file ${OUT_FILE_COMP_OUT_MERGED} \
	      --log_file ${OUT_FILE_COMP_LOG_MERGED} \
	      --min_len ${IN_MIN_LEN} \
	      --p ${IN_ALPHA} \
	      --max_calls ${IN_MAX_RECURSIVE_CALLS}
    
    ERROR=0

    for FILE in ${OUT_FILE_EVAL_UNMERGED} ${OUT_FILE_EVAL_MERGED} ${OUT_FILE_COMP_OUT_UNMERGED} ${OUT_FILE_COMP_LOG_UNMERGED} ${OUT_FILE_COMP_OUT_MERGED} ${OUT_FILE_COMP_LOG_MERGED}
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
