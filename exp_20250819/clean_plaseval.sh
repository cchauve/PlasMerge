#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
BIN_DIR=${HOME_DIR}/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output

# List of alpha values
ALPHA=$1

# Cleaning plasmid bins
logger -s "Cleaning merged bins files"
rm -f ${OUTPUT_DIR}/*/*.plaseval.merged.tsv
logger -s "Cleaning unmerged bins files"
rm -f ${OUTPUT_DIR}/*/*.plaseval.unmerged.tsv

# Cleaning eval results
logger -s "Cleaning merged eval files"
rm -f ${OUTPUT_DIR}/*/*.merged.eval.out
logger -s "Cleaning unmerged eval files"
rm -f ${OUTPUT_DIR}/*/*.unmerged.eval.out

# Cleaning comp results
for A in ${ALPHA};
do
    logger -s "Cleaning merged comp (alpha=${A}) files"
    rm -f ${OUTPUT_DIR}/*/*.merged.comp.log_${A}
    rm -f ${OUTPUT_DIR}/*/*.merged.comp.out_${A}
    logger -s "Cleaning unmerged comp (alpha=${A}) files"
    rm -f ${OUTPUT_DIR}/*/*.unmerged.comp.log_${A}
    rm -f ${OUTPUT_DIR}/*/*.unmerged.comp.out_${A}
done
