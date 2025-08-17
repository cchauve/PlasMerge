#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --account=def-chauvec
#SBATCH --job-name=plaseval_mob
#SBATCH --array=2-500
#SBATCH --output=log/plaseval/plaseval_mob_%A_%a.out
#SBATCH --error=log/plaseval/plaseval_mob_%A_%a.err

source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge/
PLASMERGE_BIN_DIR=/${HOME_DIR}/src/
PLASEVAL_BIN_DIR=/${HOME_DIR}/../PlasEval/src/
EXP_DIR=${HOME_DIR}/exp_20250802/
INPUT_PLASMERGE_DIR=${EXP_DIR}/input/
OUTPUT_PLASMERGE_DIR=${EXP_DIR}/output/
OUTPUT_PLASEVAL_DIR=${EXP_DIR}/eval_0/
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv

# Sample and input files
SAMPLE_DATA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DATA_FILE})
SAMPLE_ID=$(echo ${SAMPLE_DATA} | cut -f 1 -d ",")
ASSEMBLER=$(echo ${SAMPLE_DATA} | cut -f 2 -d ",")
EXP_ID=${SAMPLE_ID}_${ASSEMBLER}
GFA_FILE=$(echo ${SAMPLE_DATA} | cut -f 3 -d ",")
PLASCLASS_FILE=$(echo ${SAMPLE_DATA} | cut -f 5 -d ",")
PLASGRAPH_FILE=$(echo ${SAMPLE_DATA} | cut -f 6 -d ",")
MLPLASMIDS_FILE=$(echo ${SAMPLE_DATA} | cut -f 7 -d ",")
RFPLASMID_FILE=$(echo ${SAMPLE_DATA} | cut -f 8 -d ",")
MOB_FILE=$(echo ${SAMPLE_DATA} | cut -f 21 -d ",")
GT_FILE=$(echo ${SAMPLE_DATA} | cut -f 22 -d ",")

# Sample specific directories
SAMPLE_INPUT_PLASMERGE_DIR=${INPUT_PLASMERGE_DIR}/${EXP_ID}
SAMPLE_OUTPUT_PLASMERGE_DIR=${OUTPUT_PLASMERGE_DIR}/${EXP_ID}
SAMPLE_OUTPUT_PLASEVAL_DIR=${OUTPUT_PLASEVAL_DIR}/${EXP_ID}

mkdir -p ${SAMPLE_OUTPUT_PLASEVAL_DIR}

# File suffixes
RFPLASMID="rfplasmid"
PLASCLASS="plasclass"
PLASGRAPH="plasgraph2"
MLPLASMIDS="mlplasmids"
BINNING="mobrecon"
PE="plaseval"
BINNING_RFPLASMID=${BINNING}_${RFPLASMID}
BINNING_PLASCLASS=${BINNING}_${PLASCLASS}
BINNING_PLASGRAPH=${BINNING}_${PLASGRAPH}
BINNING_MLPLASMIDS=${BINNING}_${MLPLASMIDS}

# Unmerged bins files
UNMERGED_BINS_SUFFIX=tsv
UNMERGED_FILE=${SAMPLE_INPUT_PLASMERGE_DIR}/${EXP_ID}_${BINNING}.tsv
# Merged bins files
MERGED_BINS_SUFFIX=merged_bins.tsv
MERGED_RFPLASMID_FILE=${SAMPLE_OUTPUT_PLASMERGE_DIR}/${EXP_ID}_${BINNING_RFPLASMID}_${MERGED_BINS_SUFFIX}
MERGED_PLASCLASS_FILE=${SAMPLE_OUTPUT_PLASMERGE_DIR}/${EXP_ID}_${BINNING_PLASCLASS}_${MERGED_BINS_SUFFIX}
MERGED_PLASGRAPH_FILE=${SAMPLE_OUTPUT_PLASMERGE_DIR}/${EXP_ID}_${BINNING_PLASGRAPH}_${MERGED_BINS_SUFFIX}
MERGED_MLPLASMIDS_FILE=${SAMPLE_OUTPUT_PLASMERGE_DIR}/${EXP_ID}_${BINNING_MLPLASMIDS}_${MERGED_BINS_SUFFIX}

# Running Plasval
MIN_LEN=100
ALPHA=0.0
MAX_RECURSIVE_CALLS=1000000

logger -s  "#" ${EXP_ID} " running " ${PE} " on " ${BINNING} " min_len=" ${MIN_LEN} " alpha=" ${ALPHA} " max_calls=" ${MAX_RECURSIVE_CALLS}

METHOD=${BINNING}_${RFPLASMID}
logger -s  "## Evaluating " ${METHOD} "plasmid bins"
run_plaseval.sh \
    ${PLASEVAL_BIN_DIR} ${PLASMERGE_BIN_DIR} ${SAMPLE_ID} ${ASSEMBLER} ${METHOD} \
    ${GFA_FILE} ${UNMERGED_FILE} ${MERGED_RFPLASMID_FILE} \
    ${GT_FILE} \
    ${SAMPLE_OUTPUT_PLASEVAL_DIR} \
    ${MIN_LEN} ${ALPHA} ${MAX_RECURSIVE_CALLS}

METHOD=${BINNING}_${PLASCLASS}
logger -s  "## Evaluating " ${METHOD} "plasmid bins"
run_plaseval.sh \
    ${PLASEVAL_BIN_DIR} ${PLASMERGE_BIN_DIR} ${SAMPLE_ID} ${ASSEMBLER} ${METHOD} \
    ${GFA_FILE} ${UNMERGED_FILE} ${MERGED_PLASCLASS_FILE} \
    ${GT_FILE} \
    ${SAMPLE_OUTPUT_PLASEVAL_DIR} \
    ${MIN_LEN} ${ALPHA} ${MAX_RECURSIVE_CALLS}

METHOD=${BINNING}_${PLASGRAPH}
logger -s  "## Evaluating " ${METHOD} "plasmid bins"
run_plaseval.sh \
    ${PLASEVAL_BIN_DIR} ${PLASMERGE_BIN_DIR} ${SAMPLE_ID} ${ASSEMBLER} ${METHOD} \
    ${GFA_FILE} ${UNMERGED_FILE} ${MERGED_PLASGRAPH_FILE} \
    ${GT_FILE} \
    ${SAMPLE_OUTPUT_PLASEVAL_DIR} \
    ${MIN_LEN} ${ALPHA} ${MAX_RECURSIVE_CALLS}

METHOD=${BINNING}_${MLPLASMIDS}
logger -s  "## Evaluating " ${METHOD} "plasmid bins"
run_plaseval.sh \
    ${PLASEVAL_BIN_DIR} ${PLASMERGE_BIN_DIR} ${SAMPLE_ID} ${ASSEMBLER} ${METHOD} \
    ${GFA_FILE} ${UNMERGED_FILE} ${MERGED_MLPLASMIDS_FILE} \
    ${GT_FILE} \
    ${SAMPLE_OUTPUT_PLASEVAL_DIR} \
    ${MIN_LEN} ${ALPHA} ${MAX_RECURSIVE_CALLS}
