#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge/
BIN_DIR=/${HOME_DIR}/src/
EXP_DIR=${HOME_DIR}/exp_20250802/
INPUT_DIR=${EXP_DIR}/input/
OUTPUT_DIR=${EXP_DIR}/output/
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
GC_BINS_FILE=${EXP_DIR}/gc.txt

REPORT_FILE=$1
rm -f ${REPORT_FILE}
touch ${REPORT_FILE}
NB_SAMPLES=$2

RFPLASMID="rfplasmid"
PLASCLASS="plasclass"
PLASGRAPH="plasgraph2"
MLPLASMIDS="mlplasmids"
GT="ground_truth"
GP="gplascc"
MOB="mobrecon"
PBF="plasbinflow"

# File suffixes
OUTPUT_BINS_SCORES=scores.tsv
OUTPUT_MERGED_BINS=merged_bins.tsv

OK_MSG="present_non_empty"
EMPTY_MSG="empty"
MISSING_MSG="missing"

for SLURM_ARRAY_TASK_ID in $(seq 2 ${NB_SAMPLES});
do
   
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
    SAMPLE_INPUT_DIR=${INPUT_DIR}/${EXP_ID}
    SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${EXP_ID}
    
    declare -A CLASSIFICATION_INPUT_FILES=([${RFPLASMID}]=${RFPLASMID_FILE} [${PLASCLASS}]=${PLASCLASS_FILE} [${PLASGRAPH}]=${PLASGRAPH_FILE} [${MLPLASMIDS}]=${MLPLASMIDS_FILE})

    for BINNING in ${GT} ${GP} ${MOB} ${PBF};
    do
	for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
	do
	    METHOD=${BINNING}_${CLASSIFICATION}
	    CLASSIFICATION_INPUT_FILE=${CLASSIFICATION_INPUT_FILES[${CLASSIFICATION}]}
	    
	    # Converted input file name
	    FILE_UNMERGED_BINS=${SAMPLE_INPUT_DIR}/${EXP_ID}_${METHOD}.tsv
	    # PlasMerge output files
	    FILE_SCORES=${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${METHOD}_${OUTPUT_BINS_SCORES}
	    FILE_MERGED_BINS=${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${METHOD}_${OUTPUT_MERGED_BINS}

	    ERROR=0
	    REPORT_SCORES=${OK_MSG}
	    REPORT_BINS=${OK_MSG}
	    
	    if ! [ -f ${FILE_SCORES} ];
	    then
		ERROR=1
		REPORT_SCORES=${MISSING_MSG}
	    else
		FILE_SCORES_SIZE=$(stat -c%s "${FILE_SCORES}")
		if (( ${FILE_SCORES_SIZE} == 0 ));
		then
		    REPORT_SCORES=${EMPTY_MSG}
		    ERROR=1
		fi
	    fi

	    if ! [ -f ${FILE_MERGED_BINS} ];
	    then
		ERROR=1
		REPORT_BINS=${MISSING_MSG}
	    else
		FILE_MERGED_BINS_SIZE=$(stat -c%s "${FILE_MERGED_BINS}")
		if (( ${FILE_MERGED_BINS_SIZE} == 0 ));
		then
		    REPORT_BINS=${EMPTY_MSG}
		    ERROR=1
		fi
	    fi
	    
	    if (( ${ERROR} == 0 ));
	    then
		logger -s "#SUCCESS" ${EXP_ID}_${METHOD} 2>> ${REPORT_FILE}
	    else
		logger -s  "#ERROR" ${EXP_ID}_${METHOD} 2>> ${REPORT_FILE}
	    fi	  
	    logger -s "##INPUT_GFA" ${GFA_FILE} 2>> ${REPORT_FILE}
	    logger -s "##INPUT_GC" ${GC_BINS_FILE} 2>> ${REPORT_FILE}
	    logger -s "##INPUT_PLASMID_SCORES" ${CLASSIFICATION_INPUT_FILE} 2>> ${REPORT_FILE}
	    logger -s "##INPUT_UNMERGED_BINS" ${FILE_UNMERGED_BINS} 2>> ${REPORT_FILE}
	    logger -s "##OUTPUT_DIR " ${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${EXP} 2>> ${REPORT_FILE}
	    logger -s "##OUTPUT_MERGING_SCORES" ${FILE_SCORES} ${REPORT_SCORES} 2>> ${REPORT_FILE}
	    logger -s "##OUTPUT_MERGED_BINS" ${FILE_MERGED_BINS} ${REPORT_BINS} 2>> ${REPORT_FILE}
	done
    done
done

ERRORS_FILE=${REPORT_FILE}.errors
grep ERROR ${REPORT_FILE} > ${ERRORS_FILE}
HEADER=`head -1 ${DATA_FILE} `",bining,classification"
RUN_ERRORS_FILE=`echo ${DATA_FILE} | sed 's/csv/plasmerge_errors_run1.csv/g'`
echo ${HEADER} > ${RUN_ERRORS_FILE}
for BINNING in ${GT} ${GP} ${MOB} ${PBF};
do
    for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
    do
	ERRORS_FILE_METHOD=${ERRORS_FILE}.${BINNING}_${CLASSIFICATION}
	grep ${BINNING}_${CLASSIFICATION} ${ERRORS_FILE} | sed 's/_/ /g' | awk '{printf("%s,%s\n",$6,$7)}' > ${ERRORS_FILE_METHOD}
	cat ${ERRORS_FILE_METHOD} | while read line 
	do
	    SAMPLE=`grep ${line} ${DATA_FILE}`","${BINNING}","${CLASSIFICATION}
	    echo ${SAMPLE} >> ${RUN_ERRORS_FILE}
	done
	echo ${BINNING} ${CLASSIFICATION} `wc -l ${ERRORS_FILE_METHOD} | awk '{print $1}'`
	rm -f ${ERRORS_FILE_METHOD}
    done
done
rm -f ${ERRORS_FILE}
RUN_ERRORS_SAMPLES_FILE=`echo ${RUN_ERRORS_FILE} | sed 's/csv/samples.txt/g'`
cut -f1,2 -d"," ${RUN_ERRORS_FILE} | grep -v "sample" | sort -u > ${RUN_ERRORS_SAMPLES_FILE}
