# PlasMerge experiments

The experiments aim to apply PlasMerge to a set of samples for which we have
- ground truth bins,
- classification scores obtained using `plasclass,plasgraph2,mlplasmids,rfplasmid`,
- binning results obtained with these classification scores and the methods `mobrecon` (does not use classification) and `gplascc,plasbinflow`.
For each combination `ground_truth,mobrecon,gplacc,plasbinflow x plasclass,plasgraph2,mlplasmids,rfplasmid`, we run PlasMerge and then compare
the evaluation scores obtained with PlasEval for the unmerged and merged bins.


## Data files
Reusing the data file from the `20250802` experiments.
```
> cp ../exp_20250802/plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv .
> cp ../exp_20250802/plasmids_benchmarking_2025-08-02_data.txt .
```

## Python scripts
- `run_utils.py`: commands to run PlasMerge and PlasEval and to check their results;  
- `analysis_utils.py`: commands to summarize results and create figures.

## PlasMerge

Two scripts:
- `run_plasmerge.sh NB_SAMPLES`: runs PlasMerge on first `NB_SAMPLES` from `plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv`;
- results are stored in directory `output`;
- `check_plasmerge.sh NB_SAMPLES`: checks if PlasMerge results are there;
- creates files `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.<date>.[csv,txt]` (samples to rerun, log of all files).


## PlasEval

Three scripts:
- `run_plaseval_eval.sh NB_SAMPLES`: runs PlasEval in mode `eval` on first `NB_SAMPLES` from `plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv`;
- `run_plaseval_comp.sh NB_SAMPLES ALPHA`: runs PlasEval in mode `comp` on first `NB_SAMPLES` from `plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv` with `alpa=ALPHA`;
- results are stored in directory `output`, with PlasMerge results;
- `check_plaseval.sh NB_SAMPLES ALPHAS`: checks if PlasEval results are there, with `ALPHAS` being a comma-separated list of values for `alpha`.;
- creates files `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.<date>.[csv,txt]` (samples to rerun, log of all files).
