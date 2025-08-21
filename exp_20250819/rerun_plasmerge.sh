#!/bin/bash

# Data and results directories and files
HOME_DIR=/scratch/chauvec/PLASMERGE/PlasMerge
BIN_DIR=${HOME_DIR}/src
EXP_DIR=${HOME_DIR}/exp_20250819
OUTPUT_DIR=${EXP_DIR}/output
GC_BINS_FILE=${EXP_DIR}/gc.txt

DATA_FILE=$1
NB_SAMPLES=$2
RUN=$3

DATE=`date --rfc-3339=date`
mkdir -p ${OUTPUT_DIR}
SLURM_DIR=${EXP_DIR}/slurm/plasmerge/${DATE}
mkdir -p ${SLURM_DIR}
LOG_DIR=${EXP_DIR}/log/plasmerge/${DATE}
mkdir -p ${LOG_DIR}

SLURM_FILE=${SLURM_DIR}/rerun_${RUN}.sh
rm -f ${SLURM_FILE}
rm -f ${LOG_DIR}/rerun_${RUN}_*_*.out
rm -f ${LOG_DIR}/rerun_${RUN}_*_*.err
echo "#!/bin/bash" > ${SLURM_FILE}
echo "#SBATCH --time=12:00:00" >> ${SLURM_FILE}
echo "#SBATCH --mem=12G" >> ${SLURM_FILE}
echo "#SBATCH --account=def-chauvec"  >> ${SLURM_FILE}
echo "#SBATCH --job-name=plasmerge_rerun_${RUN}" >> ${SLURM_FILE}
echo "#SBATCH --array=1-${NB_SAMPLES}" >> ${SLURM_FILE}
echo "#SBATCH --output=${LOG_DIR}/rerun_${RUN}_%A_%a.out" >> ${SLURM_FILE}
echo "#SBATCH --error=${LOG_DIR}/rerun_${RUN}_%A_%a.err" >> ${SLURM_FILE}
echo "module load StdEnv/2020 gurobi/9.1.2"  >> ${SLURM_FILE}
echo "source /home/chauvec/projects/ctb-chauvec/PLASMIDS/tools/env_plasbinflow/bin/activate" >> ${SLURM_FILE}
echo "python run_utils.py plasmerge \\" >> ${SLURM_FILE}
echo "       -d ${DATA_FILE} \\" >> ${SLURM_FILE}
echo "       -i \${SLURM_ARRAY_TASK_ID} \\" >> ${SLURM_FILE}
echo "       -p ${BIN_DIR} \\" >> ${SLURM_FILE}
echo "       -gc ${GC_BINS_FILE} \\" >> ${SLURM_FILE}
echo "       -o ${OUTPUT_DIR} \\" >> ${SLURM_FILE}
echo "       -v 3.9\\" >> ${SLURM_FILE}
echo "       --rerun" >> ${SLURM_FILE}

chmod 755 ${SLURM_FILE}
sbatch ${SLURM_FILE}
