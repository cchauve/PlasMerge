#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --mem=8G
#SBATCH --account=def-chauvec
#SBATCH --job-name=plasmerge_convert
#SBATCH --array=2-1234
#SBATCH --output=log/convert/plasmerge_convert_%A_%a.out
#SBATCH --error=log/convert/plasmerge_convert_%A_%a.err

source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate
module load StdEnv/2020 gurobi/9.1.2

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge/
BIN_DIR=/${HOME_DIR}/src/
EXP_DIR=${HOME_DIR}/exp_20250802/
INPUT_DIR=${EXP_DIR}/input/
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
GC_BINS_FILE=${EXP_DIR}/gc.txt

# Sample and input files
SAMPLE_DATA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${DATA_FILE})
SAMPLE_ID=$(echo ${SAMPLE_DATA} | cut -f 1 -d ",")
ASSEMBLER=$(echo ${SAMPLE_DATA} | cut -f 2 -d ",")
EXP_ID=${SAMPLE_ID}_${ASSEMBLER}
PLASCLASS_FILE=$(echo ${SAMPLE_DATA} | cut -f 5 -d ",")
PLASGRAPH_FILE=$(echo ${SAMPLE_DATA} | cut -f 6 -d ",")
MLPLASMIDS_FILE=$(echo ${SAMPLE_DATA} | cut -f 7 -d ",")
RFPLASMID_FILE=$(echo ${SAMPLE_DATA} | cut -f 8 -d ",")
PBF_PLASCLASS_FILE=$(echo ${SAMPLE_DATA} | cut -f 13 -d ",")
PBF_PLASGRAPH_FILE=$(echo ${SAMPLE_DATA} | cut -f 14 -d ",")
PBF_MLPLASMIDS_FILE=$(echo ${SAMPLE_DATA} | cut -f 15 -d ",")
PBF_RFPLASMID_FILE=$(echo ${SAMPLE_DATA} | cut -f 16 -d ",")
GP_PLASCLASS_FILE=$(echo ${SAMPLE_DATA} | cut -f 17 -d ",")
GP_PLASGRAPH_FILE=$(echo ${SAMPLE_DATA} | cut -f 18 -d ",")
GP_MLPLASMIDS_FILE=$(echo ${SAMPLE_DATA} | cut -f 19 -d ",")
GP_RFPLASMID_FILE=$(echo ${SAMPLE_DATA} | cut -f 20 -d ",")
MOB_FILE=$(echo ${SAMPLE_DATA} | cut -f 21 -d ",")
GT_FILE=$(echo ${SAMPLE_DATA} | cut -f 22 -d ",")

# Sample specific directories
logger -s  "# 0. " ${EXP_ID} " creating directories"
SAMPLE_INPUT_DIR=${INPUT_DIR}/${EXP_ID}
rm -rf ${SAMPLE_INPUT_DIR}
mkdir ${SAMPLE_INPUT_DIR}

# File suffixes
PLASCLASS="plasclass"
PLASGRAPH="plasgraph2"
MLPLASMIDS="mlplasmids"
RFPLASMID="rfplasmid"
PBF="plasbinflow"
GP="gplascc"
MOB="mobrecon"
GT="ground_truth"
PBF_PLASCLASS=${PBF}_${PLASCLASS}
PBF_PLASGRAPH=${PBF}_${PLASGRAPH}
PBF_MLPLASMIDS=${PBF}_${MLPLASMIDS}
PBF_RFPLASMID=${PBF}_${RFPLASMID}
GP_PLASCLASS=${GP}_${PLASCLASS}
GP_PLASGRAPH=${GP}_${PLASGRAPH}
GP_MLPLASMIDS=${GP}_${MLPLASMIDS}
GP_RFPLASMID=${GP}_${RFPLASMID}

# Converted input file names
INPUT_PBF_PLASCLASS_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${PBF_PLASCLASS}.tsv
INPUT_PBF_PLASGRAPH_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${PBF_PLASGRAPH}.tsv
INPUT_PBF_MLPLASMIDS_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${PBF_MLPLASMIDS}.tsv
INPUT_PBF_RFPLASMID_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${PBF_RFPLASMID}.tsv
INPUT_GP_PLASCLASS_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${GP_PLASCLASS}.tsv
INPUT_GP_PLASGRAPH_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${GP_PLASGRAPH}.tsv
INPUT_GP_MLPLASMIDS_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${GP_MLPLASMIDS}.tsv
INPUT_GP_RFPLASMID_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${GP_RFPLASMID}.tsv
INPUT_MOB_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${MOB}.tsv
INPUT_GT_FILE=${SAMPLE_INPUT_DIR}/${EXP_ID}_${GT}.tsv

# File format conversion
logger -s  "#" ${EXP_ID} " converting input"
logger -s  "## Converting " ${PBF}+${RFPLASMID} "plasmid bins"
run_convert.sh ${BIN_DIR} ${PBF_RFPLASMID_FILE} ${INPUT_PBF_RFPLASMID_FILE} gt
logger -s  "## Converting " ${PBF}+${PLASCLASS} "plasmid bins"
run_convert.sh ${BIN_DIR} ${PBF_PLASCLASS_FILE} ${INPUT_PBF_PLASCLASS_FILE} gt
logger -s  "## Converting " ${PBF}+${PLASGRAPH} "plasmid bins"
run_convert.sh ${BIN_DIR} ${PBF_PLASGRAPH_FILE} ${INPUT_PBF_PLASGRAPH_FILE} gt
logger -s  "## Converting " ${PBF}+${MLPLASMIDS} "plasmid bins"
run_convert.sh ${BIN_DIR} ${PBF_MLPLASMIDS_FILE} ${INPUT_PBF_MLPLASMIDS_FILE} gt
logger -s  "## Converting " ${GP}+${RFPLASMID} "plasmid bins"
run_convert.sh ${BIN_DIR} ${GP_RFPLASMID_FILE} ${INPUT_GP_RFPLASMID_FILE} gt
logger -s  "## Converting " ${GP}+${PLASCLASS} "plasmid bins"
run_convert.sh ${BIN_DIR} ${GP_PLASCLASS_FILE} ${INPUT_GP_PLASCLASS_FILE} gt
logger -s  "## Converting " ${GP}+${PLASGRAPH} "plasmid bins"
run_convert.sh ${BIN_DIR} ${GP_PLASGRAPH_FILE} ${INPUT_GP_PLASGRAPH_FILE} gt
logger -s  "## Converting " ${GP}+${MLPLASMIDS} "plasmid bins"
run_convert.sh ${BIN_DIR} ${GP_MLPLASMIDS_FILE} ${INPUT_GP_MLPLASMIDS_FILE} gt
logger -s  "## Converting " ${MOB} "plasmid bins"
run_convert.sh ${BIN_DIR} ${MOB_FILE} ${INPUT_MOB_FILE} gt
logger -s  "## Converting " ${GT} "plasmid bins"
run_convert.sh ${BIN_DIR} ${GT_FILE} ${INPUT_GT_FILE} gt
