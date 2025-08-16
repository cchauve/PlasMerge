#!/bin/bash

source ../../env_plaseval/bin/activate

logger -s "Creating CSV file with alpha=0.5"
python PlasEval_utils.py csv \
       plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv \
       eval \
       analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv \
       -n 499 \
       -v \
       > analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.NA.txt

logger -s "Creating CSV file with alpha=0"
python PlasEval_utils.py csv \
       plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv \
       eval_0 \
       analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results_0.csv \
       -n 499 \
       -v \
       > analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results_0.NA.txt

logger -s "Creating plots with alpha=0.5"
./create_scatter_plots.sh analysis/figures analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv
./create_differemce_plots.sh analysis/figures analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv

logger -s "Creating plots with alpha=0"
./create_scatter_plots.sh analysis/figures_0 analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results_0.csv
./create_difference_plots.sh analysis/figures_0 analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results_0.csv
