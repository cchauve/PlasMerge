#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
BIN_DIR=${HOME_DIR}/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output
DATA_FILE=${EXP_DIR}/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
GC_BINS_FILE=${EXP_DIR}/gc.txt

NB_SAMPLES=$1

DATE=`date --rfc-3339=date`
mkdir -p ${OUTPUT_DIR}
SLURM_DIR=${EXP_DIR}/slurm/plasmerge/${DATE}
mkdir -p ${SLURM_DIR}
LOG_DIR=${EXP_DIR}/log/plasmerge/${DATE}
mkdir -p ${LOG_DIR}

for BINNING in "ground_truth" "gplascc" "mobrecon" "plasbinflow";
do
    for CLASSIFICATION in "plasclass" "plasgraph2" "rfplasmid" "mlplasmids";
    do
	SLURM_FILE=${SLURM_DIR}/run_${BINNING}_${CLASSIFICATION}.sh
	rm -f ${SLURM_FILE}
	rm -f ${LOG_DIR}/${BINNING}_${CLASSIFICATION}_*_*.out
	rm -f ${LOG_DIR}/${BINNING}_${CLASSIFICATION}_*_*.err
	echo "#!/bin/bash" > ${SLURM_FILE}
	echo "#SBATCH --time=8:00:00" >> ${SLURM_FILE}
	echo "#SBATCH --mem=8G" >> ${SLURM_FILE}
	echo "#SBATCH --account=def-chauvec"  >> ${SLURM_FILE}
	echo "#SBATCH --job-name=plasmerge_${BINNING}_${CLASSIFICATION}" >> ${SLURM_FILE}
	echo "#SBATCH --array=1-${NB_SAMPLES}" >> ${SLURM_FILE}
	echo "#SBATCH --output=${LOG_DIR}/${BINNING}_${CLASSIFICATION}_%A_%a.out" >> ${SLURM_FILE}
	echo "#SBATCH --error=${LOG_DIR}/${BINNING}_${CLASSIFICATION}_%A_%a.err" >> ${SLURM_FILE}
	echo "module load StdEnv/2020 gurobi/9.1.2"  >> ${SLURM_FILE}
	echo "source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate" >> ${SLURM_FILE}
	echo "python run_utils.py plasmerge \\" >> ${SLURM_FILE}
	echo "       -d ${DATA_FILE} \\" >> ${SLURM_FILE}
	echo "       -i \${SLURM_ARRAY_TASK_ID} \\" >> ${SLURM_FILE}
	echo "       -b ${BINNING} \\" >> ${SLURM_FILE}
	echo "       -c ${CLASSIFICATION} \\" >> ${SLURM_FILE}
	echo "       -p ${BIN_DIR} \\" >> ${SLURM_FILE}
	echo "       -gc ${GC_BINS_FILE} \\" >> ${SLURM_FILE}
	echo "       -o ${OUTPUT_DIR} \\" >> ${SLURM_FILE}
	echo "       -v 3.9" >> ${SLURM_FILE}

	chmod 755 ${SLURM_FILE}
	sbatch ${SLURM_FILE}
    done
done
