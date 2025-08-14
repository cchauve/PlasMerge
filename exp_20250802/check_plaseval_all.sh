#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge/
BIN_DIR=/${HOME_DIR}/src/
EXP_DIR=${HOME_DIR}/exp_20250802/
#OUTPUT_DIR=${EXP_DIR}/eval/
OUTPUT_DIR=${EXP_DIR}/eval_0/
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv

REPORT_FILE=$1
rm -f ${REPORT_FILE}
touch ${REPORT_FILE}
NB_SAMPLES=$2

# Methods
RFPLASMID="rfplasmid"
PLASCLASS="plasclass"
PLASGRAPH="plasgraph2"
MLPLASMIDS="mlplasmids"
GT="ground_truth"
GP="gplascc"
MOB="mobrecon"
PBF="plasbinflow"
MERGED="merged"
UNMERGED="unmerged"

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
    SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${EXP_ID}
    
    for BINNING in ${GT} ${GP} ${MOB} ${PBF};
    do
	for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
	do
	    for BINS in ${UNMERGED} ${MERGED};
	    do
		METHOD=${BINNING}_${CLASSIFICATION}
		INPUT_BINS_FILE=${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${METHOD}.${BINS}.tsv
		FILE_EVAL=${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${METHOD}.${BINS}.eval.out
		FILE_COMP_OUT=${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${METHOD}.${BINS}.comp.out
		FILE_COMP_LOG=${SAMPLE_OUTPUT_DIR}/${EXP_ID}_${METHOD}.${BINS}.comp.log
	    
		ERROR=0
		REPORT_INPUT=${OK_MSG}
		REPORT_EVAL=${OK_MSG}
		REPORT_COMP_OUT=${OK_MSG}
		REPORT_COMP_LOG=${OK_MSG}
	    
		if ! [ -f ${INPUT_BINS_FILE} ];
		then
		    ERROR=1
		    REPORT_INPUT=${MISSING_MSG}
		else
		    FILE_INPUT_SIZE=$(stat -c%s "${INPUT_BINS_FILE}")
		    if (( ${FILE_INPUT_SIZE} == 0 ));
		    then
			REPORT_INPUT=${EMPTY_MSG}
			ERROR=1
		    fi
		fi

		if ! [ -f ${FILE_EVAL} ];
		then
		    ERROR=1
		    REPORT_EVAL=${MISSING_MSG}
		else
		    FILE_EVAL_SIZE=$(stat -c%s "${FILE_EVAL}")
		    if (( ${FILE_EVAL_SIZE} == 0 ));
		    then
			REPORT_EVAL=${EMPTY_MSG}
			ERROR=1
		    fi
		fi

		if ! [ -f ${FILE_COMP_OUT} ];
		then
		    ERROR=1
		    REPORT_COMP_OUT=${MISSING_MSG}
		else
		    FILE_COMP_OUT_SIZE=$(stat -c%s "${FILE_COMP_OUT}")
		    if (( ${FILE_COMP_OUT_SIZE} == 0 ));
		    then
			REPORT_COMP_OUT=${EMPTY_MSG}
			ERROR=1
		    fi
		fi

		if ! [ -f ${FILE_COMP_LOG} ];
		then
		    ERROR=1
		    REPORT_COMP_LOG=${MISSING_MSG}
		else
		    FILE_COMP_LOG_SIZE=$(stat -c%s "${FILE_COMP_LOG}")
		    if (( ${FILE_COMP_LOG_SIZE} == 0 ));
		    then
			REPORT_COMP_LOG=${EMPTY_MSG}
			ERROR=1
		    fi
		fi
	    
		if (( ${ERROR} == 0 ));
		then
		    logger -s "#SUCCESS" ${EXP_ID}_${METHOD}_${BINS} 2>> ${REPORT_FILE}
		else
		    logger -s  "#ERROR" ${EXP_ID}_${METHOD}_${BINS} 2>> ${REPORT_FILE}
		fi	  
		logger -s "##INPUT_BINS" ${INPUT_BINS_FILE} ${REPORT_INPUT} 2>> ${REPORT_FILE}
		logger -s "##INPUT_GROUND_TRUTH" ${GT_FILE} 2>> ${REPORT_FILE}
		logger -s "##OUTPUT_DIR " ${SAMPLE_OUTPUT_DIR} 2>> ${REPORT_FILE}
		logger -s "##OUTPUT_EVAL" ${FILE_EVAL} ${REPORT_EVAL} 2>> ${REPORT_FILE}
		logger -s "##OUTPUT_COMP_OUT" ${FILE_COMP_OUT} ${REPORT_COMP_OUT} 2>> ${REPORT_FILE}
		logger -s "##OUTPUT_COMP_LOG" ${FILE_COMP_LOG} ${REPORT_COMP_LOG} 2>> ${REPORT_FILE}
	    done
	done
    done
done

# ERRORS_FILE=${REPORT_FILE}.errors
# grep ERROR ${REPORT_FILE} > ${ERRORS_FILE}
# HEADER=`head -1 ${DATA_FILE} `",bining,classification,bins"
# RUN_ERRORS_FILE=`echo ${DATA_FILE} | sed 's/csv/plaseval_errors_run1.csv/g'`
# echo ${HEADER} > ${RUN_ERRORS_FILE}
# for BINNING in ${GT} ${GP} ${MOB} ${PBF};
# do
#     for CLASSIFICATION in ${RFPLASMID} ${PLASCLASS} ${PLASGRAPH} ${MLPLASMIDS};
#     do
#         for BINS in ${UNMERGED} ${MERGED};
# 	do
# 	    ERRORS_FILE_METHOD=${ERRORS_FILE}.${BINNING}_${CLASSIFICATION}_${BINS}
#             grep ${BINNING}_${CLASSIFICATION}_${BINS} ${ERRORS_FILE} | sed 's/_/ /g' | awk '{printf("%s,%s\n",$6,$7)}' > ${ERRORS_FILE_METHOD}
#             cat ${ERRORS_FILE_METHOD} | while read line
#             do
# 		SAMPLE=`grep ${line} ${DATA_FILE}`","${BINNING}","${CLASSIFICATION}","${BINS}
# 		echo ${SAMPLE} >> ${RUN_ERRORS_FILE}
#             done
#             echo ${BINNING} ${CLASSIFICATION} ${BINS} `wc -l ${ERRORS_FILE_METHOD} | awk '{print $1}'`
#             rm -f ${ERRORS_FILE_METHOD}
# 	done
#     done
# done
# rm -f ${ERRORS_FILE}
# RUN_ERRORS_SAMPLES_FILE=`echo ${RUN_ERRORS_FILE} | sed 's/csv/samples.txt/g'`
# cut -f1,2 -d"," ${RUN_ERRORS_FILE} | grep -v "sample" | sort -u > ${RUN_ERRORS_SAMPLES_FILE}
