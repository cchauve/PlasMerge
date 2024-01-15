#!/bin/bash

source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate
module load gurobi

SAMPLES='SAMN32247302 SAMN32247345 SAMN32247425 SAMN32247519 SAMN32247522'
SOURCES='gt mob gp pbf'

for SAMPLE in ${SAMPLES}
do
    for SOURCE in ${SOURCES}
    do
	python src/plasmerge.py \
	       -s ${SAMPLE} \
	       -a test/gfas/${SAMPLE}.assembly.gfa.gz \
	       -p test/scores/${SAMPLE}.scores.tsv \
	       -b test/pls_bins/${SAMPLE}.${SOURCE}.tsv \
	       -r ${SOURCE} \
	       -g test/gc_intervals.txt \
	       -d test/results/model \
	       -f test/results/${SOURCE}/${SAMPLE}.${SOURCE}.tsv \
	       > test/results/${SAMPLE}.${SOURCE}.plasmerge.out
    done
done
