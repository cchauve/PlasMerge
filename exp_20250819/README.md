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
drwxr-sr-x 2 chauvec def-chauvec  131072 Sep  7 18:00 figures/
-rw-r--r-- 1 chauvec def-chauvec   56740 Sep  7 18:06 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.0.2025-09-07.NA.txt
-rw-r--r-- 1 chauvec def-chauvec 2924049 Sep  7 18:06 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.0.2025-09-07.csv
-rw-r--r-- 1 chauvec def-chauvec   56740 Sep  7 18:05 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.5.2025-09-07.NA.txt
-rw-r--r-- 1 chauvec def-chauvec 3408589 Sep  7 18:05 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.0.5.2025-09-07.csv
-rw-r--r-- 1 chauvec def-chauvec   56740 Sep  7 18:06 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.1.0.2025-09-07.NA.txt
-rw-r--r-- 1 chauvec def-chauvec 3100317 Sep  7 18:06 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.1.0.2025-09-07.csv
```

## 2025-09-07: output files nomenclature

### Sample-specific files
All the files generated for a pair sample+assembler `<SAMPLE_ASSEMBLER>`are in the directory `output/<SAMPLE_ASSEMBLER>`.
These files (if all experiments succeeded, othewise some file can be missing) are:
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.plasmerge.unmerged.tsv`: unmerged bins obtained by running `<BINNING_METHOD>` with plasmd scores obtained with `<CLASSIF_METHOD>`, in PlasMerge format;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.plaseval.unmerged.tsv`: same bins than file above, but in PlasEval format;  
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.plasmerge.merged.tsv`: merged bins (by PlasMerge) from the unmerged bins file above, in PlasMerge format;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.plaseval.merged.tsv`: same bins than file above, but in PlasEval format;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.plasmerge.merging_scores.tsv`: merging scores generated by PlasMerge;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.unmerged.eval.out`: result of PlasEval in `eval` mode on unmerged bins, compared to gound truth bins;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.merged.eval.out`: result of PlasEval in `eval` mode on merged bins, compared to gound truth bins;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.unmerged.comp.out_<ALPHA>`: result of PlasEval in `comp` mode with `alpha` set to `<ALPHA>` on unmerged bins, compared to gound truth bins;
- `<SAMPLE_ASSEMBLER>.<BINNING_METHOD>_<CLASSIF_METHOD>.merged.comp.out_<ALPHA>`: result of PlasEval in `comp` mode with `alpha` set to `<ALPHA>` on merged bins, compared to gound truth bins.

The plasmid binning methods are `ground_truth` (true bins), `gplascc`, `plasbinflow` and `mobrecon`.
The classification methods are `plasclass`, `plasgraph2`, `mlplasmids`, `rfplasmid`.
Three alpha values were considered: `0.0` (count cuts and joins only), `0.5` and `1.0`.

### PlasEval summary files
For every value of `alpha`, all statistcs on
- number of bins (`nb_bins`), of contigs (`nb_ctgs`), length of contigs (`len_ctgs`) in bins;
- PlasEval `eval` mode (`Precision`, `Recall`, `F1`: do not depend on `alpha`), in weighted (`w`) and unweighted (`u`) mode;
- PlasEval `comp` mode (`Dissimilarity`, `Cuts`, `Joins`, `Extra_ctgs`, `Missing_ctgs`) in normalized (`n`) and unnormalized (`u`) mode;  
are recorded in the file `analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval.<ALPHA>.<DATE>.tsv`.

### Figures 
For each statistic listed in the previous paragraph, we generate figures allowing to compare the results between the unmerged bins and the merged bins; these comparisons are shown in two forms, a scatter plot and a violin plot of the distribution of the differences of values (`merged` - `unmerged`), both over each combination `<BINNING_METHOD>,<CLASSIF_METHOD>` (one scatter plot file and one violin plot file per combination) and over all combinations aggregated (again one scatter plot file and one difference file).

For bins and contigs statistics `nb_bins,nb_ctgs,len_ctgs`, and a combination `<BINNING_METHOD>,<CLASSIF_METHOD>`, the scatter plot and violin plot files for a given statistic `<STAT>` are
- `analysis/figures/[scatter,difference]_<STAT>_<BINNING METHOD>_<CLASSIF_METHOD>.png`.  
The files for statistics aggregated over all methods combinations are
- `analysis/figures/[scatter,difference]_<STAT>_aggregated.png`.  

For a statistic in `Precision,Recall,F1` and a mode `<MODE>` (weighted `w`/unweighted `u`) and a combination `<BINNING_METHOD>,<CLASSIF_METHOD>`, the scatter plot and violin plot files for the given statistic `<STAT>` are
- `analysis/figures/[scatter,difference]_<MODE>.<STAT>_<BINNING METHOD>_<CLASSIF_METHOD>.png`.  
The files for statistics aggregated over all methods combinations are
- `analysis/figures/[scatter,difference]_<MODE>.<STAT>_aggregated.png`.  

For a statistic in `Dissimilarity,Cuts,Joins,Extra_ctgs,Missing_ctgs` and a mode `<MODE>` (normalized `n`/unnormalized `u`), obtained with `alpha` value `<ALPHA>`, and a combination `<BINNING_METHOD>,<CLASSIF_METHOD>`, the scatter plot and violin plot files for the given statistic `<STAT>` are
- `analysis/figures/[scatter,difference]_<MODE>.<STAT>_<ALPHA>_<BINNING METHOD>_<CLASSIF_METHOD>.png`.  
The files for statistics aggregated over all methods combinations are
- `analysis/figures/[scatter,difference]_<MODE>.<STAT>_<ALPHA>_aggregated.png`.  

Note there are also figures for some joint statstics: `Extra_ctgs,Missing_ctgs` and `Cuts,Joins`.
