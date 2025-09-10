#!/bin/bash

source /scratch/chauvec/PLASMERGE/env_plaseval/bin/activate

DATA_FILE=plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
RES_FILE_05=analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.5.2025-09-09.csv
NA_FILE_05=analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.5.2025-09-09.NA.txt
RES_FILE_0=analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.0.2025-09-09.csv
NA_FILE_0=analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.0.2025-09-09.NA.txt
RES_FILE_1=analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.1.0.2025-09-09.csv
NA_FILE_1=analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.1.0.2025-09-09.NA.txt


# Creating results files
logger -s "Creating CSV file with alpha=0.5"
python analysis_utils.py csv ${DATA_FILE} output/ ${RES_FILE_05} -n 500 -a 0.5 > ${NA_FILE_05}
logger -s "Creating CSV file with alpha=0.0"
python analysis_utils.py csv ${DATA_FILE} output/ ${RES_FILE_0}  -n 500 -a 0.0 > ${NA_FILE_0}
logger -s "Creating CSV file with alpha=1.0"
python analysis_utils.py csv ${DATA_FILE} output/ ${RES_FILE_1}  -n 500 -a 1.0 > ${NA_FILE_1}

# Bins stats plots
logger -s "Bins statistics, scatter"
./create_plots.sh bins scatter    analysis/figures ${RES_FILE_05} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" NA NA
logger -s "Bins statistics, difference"
./create_plots.sh bins difference analysis/figures ${RES_FILE_05} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" NA NA

# Eval plots
logger -s "Eval, scatter"
./create_plots.sh eval scatter    analysis/figures ${RES_FILE_05} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" NA "u w"
logger -s "Eval, difference"
./create_plots.sh eval difference analysis/figures ${RES_FILE_05} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" NA "u w"

# Comp plots, alpha=0.5
logger -s "Comp, alpha=0.5, scatter"
./create_plots.sh comp scatter    analysis/figures ${RES_FILE_05} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" 0.5 "u n"
logger -s "Comp, alpha=0.5, difference"
./create_plots.sh comp difference analysis/figures ${RES_FILE_05} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" 0.5 "u n"

# Comp plots, alpha=1.0
logger -s "Comp, alpha=1.0, scatter"
./create_plots.sh comp scatter    analysis/figures ${RES_FILE_1} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" 1.0 "u n"
logger -s "Comp, alpha=1.0, difference"
./create_plots.sh comp difference analysis/figures ${RES_FILE_1} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" 1.0 "u n"

# Comp plots, alpha=0.0
logger -s "Comp, alpha=0.0, scatter"
./create_plots.sh comp scatter    analysis/figures ${RES_FILE_0} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" 0.0 "u n"
logger -s "Comp, alpha=0.0, difference"
./create_plots.sh comp difference analysis/figures ${RES_FILE_0} "ground_truth gplascc plasbinflow mobrecon" "plasclass plasgraph2 mlplasmids rfplasmid" 0.0 "u n"

