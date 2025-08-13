#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --mem=8G
#SBATCH --account=def-chauvec
#SBATCH --job-name=plasmerge_check_data
#SBATCH --output=log/plasmerge_check_data.out
#SBATCH --error=log/plasmerge_check_data.err

# Data files
DATA_FILE=plasmids_benchmarking_2025-08-02_data.csv
DATA_DESC_FILE=plasmids_benchmarking_2025-08-02_data.txt

# Parameters of data
## number of samples: first line of file is the header
OFFSET_SAMPLES=1
NB_SAMPLES=$(( `wc -l < ${DATA_FILE}` - ${OFFSET_SAMPLES} ))
## first 2 fields are sample information, not files
FIRST_FILE_FIELD=3
## number of fields in data file: no header
OFFSET_FIELDS_DESC=0
NB_FIELDS=$(( `wc -l < ${DATA_DESC_FILE}` - ${OFFSET_FIELDS_DESC} ))

# Output
ERRORS_FILE=plasmids_benchmarking_2025-08-02_data.errors.txt
rm -f ${ERRORS_FILE}
touch ${ERRORS_FILE}
OUT_FILE=plasmids_benchmarking_2025-08-02_data.filtered.csv
rm -f ${OUT_FILE}
touch ${OUT_FILE}
head -1 ${DATA_FILE} > ${OUT_FILE}

# Loop over samples
for i in $(seq 1 ${NB_SAMPLES});
do
    DATA_OK=1
    # Sample information
    SAMPLE_LINE=$(( ${i} + ${OFFSET_SAMPLES} ))
    SAMPLE_DATA=$(sed -n "${SAMPLE_LINE}p" ${DATA_FILE})
    SAMPLE_ID=$(echo ${SAMPLE_DATA} | cut -f 1 -d ",")
    ASSEMBLER=$(echo ${SAMPLE_DATA} | cut -f 2 -d ",")
    EXP_ID=${SAMPLE_ID}_${ASSEMBLER}
    
    # Loop over plasmerge-relevant file fields for a given sample; first 2 fields are sample information
    # $(seq ${FIRST_FILE_FIELD} ${NB_FIELDS});
    for j in 3 4 5 6 7 8 13 14 15 16 17 18 19 20 21 22;
    do
	FILE=$(echo ${SAMPLE_DATA} | cut -f ${j} -d ",")
	# File description
	k=$(( ${j} + ${OFFSET_FIELDS_DESC} ))
	FILE_DESC=$(sed -n "${k}p" ${DATA_DESC_FILE} | cut -f 2 -d " ")
	MSG="${EXP_ID} ${FILE_DESC}"
	if [ -z "${FILE}" ];
	then
	    echo "WARNING.EMPTY_FIELD" ${MSG} >> ${ERRORS_FILE}
	    DATA_OK=0
	else
	    if ! [ -f "${FILE}" ];
	    then
		echo "ERROR.MISSING_FILE" ${MSG} ${FILE} >> ${ERRORS_FILE}
		DATA_OK=0
	    else
		FILE_SIZE=$(stat -c%s "${FILE}")
		if (( ${FILE_SIZE} == 0 ));
		then
		    echo "ERROR.EMPTY_FILE" ${MSG} ${FILE} >> ${ERRORS_FILE}
		    DATA_OK=0
		fi
	    fi
	fi
    done
    if [[ "${DATA_OK}" == 1 ]]; then
	echo ${SAMPLE_DATA} >> ${OUT_FILE}
    fi
done
