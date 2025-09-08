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

## 2025-08-21
Running PlasMerge.
```
> ./run_plasmerge.sh 500
```

## 2025-08-22
Checking PlasMerge results.
```
> ./check_plasmerge.sh 500 1
ground_truth plasclass 2
ground_truth plasgraph2 3
ground_truth rfplasmid 2
ground_truth mlplasmids 2
gplascc plasclass 64
gplascc plasgraph2 36
gplascc rfplasmid 5
gplascc mlplasmids 58
mobrecon plasclass 4
mobrecon plasgraph2 1
mobrecon rfplasmid 1
mobrecon mlplasmids 1
plasbinflow plasclass 9
plasbinflow plasgraph2 4
plasbinflow rfplasmid 1
plasbinflow mlplasmids 1
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.1.2025-08-22.csv
194 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.1.2025-08-22.csv
```
Only 194 experiments did not finish properly, most of them due `gplascc`.

## 2025-08-23
Running PlasEval in `eval` mode (note: `comp` mode has a bug) for unmerged bins, then, once done, merged bins.
```
> ./run_plaseval_eval.sh 500 unmerged
> ./run_plaseval_eval.sh 500 merged
```

## 2025-09-05
Running PlasEval in `comp` mode (note: `comp` mode bug fixed) for unmerged bins, then, once done, merged bins, with `alpha=05`, then with `alpha=1`.
```
> ./run_plaseval_comp.sh 500 0.5 unmerged
> ./run_plaseval_comp.sh 500 0.5 merged
> ./run_plaseval_comp.sh 500 1 unmerged
> ./run_plaseval_comp.sh 500 1 merged
```

## 2025-09-07
Checking PlasEval in `eval` mode.
```
> ./check_plaseval.sh 500 eval 0 1
ground_truth plasclass unmerged 0
...
plasbinflow mlplasmids merged 1
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.eval.0.1.2025-09-07.csv
194 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.eval.0.1.2025-09-07.csv
```
No PlasEval `eval` experiment did fail, all missig results are due to missing PlasMerge results.  

Checking PlasEval in `comp` mode with `alpha=0.5`.

```
> ./check_plaseval.sh 500 comp 0.5 1
ground_truth plasclass unmerged 0
...
plasbinflow mlplasmids merged 1
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.5.1.2025-09-07.csv
194 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.5.1.2025-09-07.csv
```
Same as above. Checking now with `alpha=1.0`.
```
> ./check_plaseval.sh 500 comp 1.0 1
ground_truth plasclass unmerged 0
...
plasbinflow mlplasmids merged 1
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.1.0.1.2025-09-07.csv
194 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.1.0.1.2025-09-07.csv
```

## 2025-09-07
Correcting mistake: it should have been with `alpha=0.0` and not `alpha=1.0`.
```
> ./run_plaseval_comp.sh 500 0.0 unmerged
> ./run_plaseval_comp.sh 500 0.0 merged
```
Checking PlasEval in mode `comp` with `alpha=0.0`.
```
> ./check_plaseval.sh 500 comp 0.0 1
ground_truth plasclass unmerged 0
...
plasbinflow mlplasmids merged 1
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.0.1.2025-09-07.csv
194 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.0.1.2025-09-07.csv
```
Same as above.

## 2025-09-07
Creating `csv` results files and figures.
```
> ll plasmids_benchmarking*
-rw-r--r-- 1 chauvec def-chauvec  3224700 Aug 19 10:59 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
-rw-r--r-- 1 chauvec def-chauvec   514553 Sep  7 16:58 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.0.1.2025-09-07.csv
-rw-r--r-- 1 chauvec def-chauvec 19386350 Sep  7 16:58 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.0.1.2025-09-07.txt
-rw-r--r-- 1 chauvec def-chauvec   514553 Sep  7 09:18 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.5.1.2025-09-07.csv
-rw-r--r-- 1 chauvec def-chauvec 19386350 Sep  7 09:18 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.0.5.1.2025-09-07.txt
-rw-r--r-- 1 chauvec def-chauvec   514553 Sep  7 09:21 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.1.0.1.2025-09-07.csv
-rw-r--r-- 1 chauvec def-chauvec 19386350 Sep  7 09:21 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.comp.1.0.1.2025-09-07.txt
-rw-r--r-- 1 chauvec def-chauvec   514553 Sep  7 09:14 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.eval.0.1.2025-09-07.csv
-rw-r--r-- 1 chauvec def-chauvec 15937882 Sep  7 09:14 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.eval.0.1.2025-09-07.txt
-rw-r--r-- 1 chauvec def-chauvec   512018 Aug 22 08:50 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.1.2025-08-22.csv
-rw-r--r-- 1 chauvec def-chauvec  9590218 Aug 22 08:41 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge.1.2025-08-22.txt
-rw-r--r-- 1 chauvec def-chauvec      345 Aug 19 10:59 plasmids_benchmarking_2025-08-02_data.txt
> nohup ./create_plots_all.sh &
> grep -v flexiblas nohup.out > log/create_plots_all.log
> rm nohup.out
> ll analysis/

```
