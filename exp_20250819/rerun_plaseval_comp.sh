#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
PLASMERGE_BIN_DIR=${HOME_DIR}/src
PLASEVAL_BIN_DIR=${HOME_DIR}/../PlasEval/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output

DATA_FILE=$1
NB_SAMPLES=$2
ALPHA=$3
RUN=$4
MIN_LEN=100
MAX_CALLS=1000000

DATE=`date --rfc-3339=date`
mkdir -p ${OUTPUT_DIR}
SLURM_DIR=${EXP_DIR}/slurm/plaseval/${DATE}
mkdir -p ${SLURM_DIR}
LOG_DIR=${EXP_DIR}/log/plaseval/${DATE}
mkdir -p ${LOG_DIR}

PREFIX="comp_rerun"

SLURM_FILE=${SLURM_DIR}/${PREFIX}_${RUN}_${ALPHA}.sh
rm -f ${SLURM_FILE}
rm -f ${LOG_DIR}/${PREFIX}_${RUN}_${ALPHA}_*_*.out
rm -f ${LOG_DIR}/${PREFIX}_${RUN}_${ALPHA}_*_*.err
echo "#!/bin/bash" > ${SLURM_FILE}
echo "#SBATCH --time=12:00:00" >> ${SLURM_FILE}
echo "#SBATCH --mem=8G" >> ${SLURM_FILE}
echo "#SBATCH --account=def-chauvec"  >> ${SLURM_FILE}
echo "#SBATCH --job-name=plaseval_${PREFIX}_${RUN}" >> ${SLURM_FILE}
echo "#SBATCH --array=1-${NB_SAMPLES}" >> ${SLURM_FILE}
echo "#SBATCH --output=${LOG_DIR}/${PREFIX}_${RUN}_%A_%a.out" >> ${SLURM_FILE}
echo "#SBATCH --error=${LOG_DIR}/${PREFIX}_${RUN}_%A_%a.err" >> ${SLURM_FILE}
echo "source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate" >> ${SLURM_FILE}
echo "python run_utils.py plaseval \\" >> ${SLURM_FILE}
echo "       -d ${DATA_FILE} \\" >> ${SLURM_FILE}
echo "       -i \${SLURM_ARRAY_TASK_ID} \\" >> ${SLURM_FILE}
echo "       -pm ${PLASMERGE_BIN_DIR} \\" >> ${SLURM_FILE}
echo "       -pe ${PLASEVAL_BIN_DIR} \\" >> ${SLURM_FILE}
echo "       -a ${ALPHA} \\" >> ${SLURM_FILE}
echo "       -m comp \\" >> ${SLURM_FILE}
echo "       -ml ${MIN_LEN} \\" >> ${SLURM_FILE}
echo "       -mc ${MAX_CALLS} \\" >> ${SLURM_FILE}
echo "       -o ${OUTPUT_DIR} \\" >> ${SLURM_FILE}
echo "       -v 3.9 \\" >> ${SLURM_FILE}
echo "       --rerun" >> ${SLURM_FILE}

chmod 755 ${SLURM_FILE}
sbatch ${SLURM_FILE}
