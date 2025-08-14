# PlasMerge experiments

The experiments aim to apply PlasMerge to a set of samples for which we have  
- ground truth bins,  
- classification scores obtained using `plasclass,plasgraph2,mlplasmids,rfplasmid`,  
- binning results obtained with these classification scores and the methods `mobrecon` (does not use classification) and `gplascc,plasbinflow`.  
For each combination `ground_truth,mobrecon,gplacc,plasbinflow x plasclass,plasgraph2,mlplasmids,rfplasmid`, we run PlasMerge and then compare
the evaluation scores obtained with PlasEval for the unmerged and merged bins.


## Data files
Downloading data file from Aniket.
```
> wget https://github.com/acme92/benchmarking/blob/benchmarking_paths/plasmerge_experiments_paths.csv
> mv plasmerge_experiments_paths.csv plasmids_benchmarking_2025-08-02_data.csv
> head -1 plasmids_benchmarking_2025-08-02_data.csv \
  | sed 's/,/ /g' \
  | awk '{for(i=1;i<=NF;i++){printf("%d. %s\n",i,$i)}}' \
  > plasmids_benchmarking_2025-08-02_data.txt
```

## Checking data
Checking data to filter out samples for which some file is missing.
```
> sbatch check_data.sh
```
Detected errors reported in `plasmids_benchmarking_2025-08-02_data.errors.txt`.
High-level statistics on whic kinds of files are missing.
```
> wc -l plasmids_benchmarking_2025-08-02_data.csv
2483
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.csv
1234
> cut -f 2 -d " " plasmids_benchmarking_2025-08-02_data.errors.txt \
  | grep SKESA \
  | sort -u \
  | wc -l
675
> cut -f 2 -d " " plasmids_benchmarking_2025-08-02_data.errors.txt \
  | grep UNICYCLER \
  | sort -u \
  | wc -l
574
> cut -f 3 -d " " plasmids_benchmarking_2025-08-02_data.errors.txt \
  | sort \
  | awk \
  'BEGIN{TXT="";NB=0} \
  {if($1!=TXT){printf("%s %d\n",TXT,NB);TXT=$1;} else{NB++;}}\
  END{printf("%s %d\n",TXT,NB);}'
gp_mlpl_bins 1107
gp_plcl_bins 1317
gp_plgr_bins 1779
gp_rfpl_bins 2116
input_pbf_mlpl 3071
mob_bins 3072
pbf_mlpl_bins 4028
pbf_plcl_bins 4035
pbf_plgr_bins 4047
pbf_rfpl_bins 4055
```
For half of the samples, at least one file is missing.  
Both assemblers are involved in errors.  
All methods sow many errors, with `mlplasmids` failing a very large number of times.   

## Randomizing input
We will randomize the error-free samples an process 499 of them. 
```
> head -1 plasmids_benchmarking_2025-08-02_data.filtered.csv \
  > plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
> grep -v species_sample plasmids_benchmarking_2025-08-02_data.filtered.csv \
  | shuf >> plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
1234 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
```
From now on the data file used is `plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv`.  
We will process the first 499 samples.
```
> head -500 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv | grep -c abau
92
> head -500 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv | grep -c ecol
281
> head -500 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv | grep -c efae
47
> head -500 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv | grep -c kpne
79
```
High proportion of *E. coli* samples.


## PlasMerge
Converting the data into PlasMerge format.
```
> sbatch run_convert_all.sh
```
No error.

Processing the first 499 randomized samples.
```
> sbatch run_plasmerge_gt_all.sh
> sbatch run_plasmerge_gp_all.sh
> sbatch run_plasmerge_pbf_all.sh
. sbatch run_plasmerge_mob_all.sh
```
Redone on `2025-08-09` due to inconsistency in output directories between the scripts.  
Finished on `2028-08-11`.  
Checking results and recording success/errors.    
```
> ./check_plasmerge_all.sh report_plasmerge_run1_20250811.txt 500
ground_truth rfplasmid 2
ground_truth plasclass 2
ground_truth plasgraph2 2
ground_truth mlplasmids 3
gplascc rfplasmid 3
gplascc plasclass 35
gplascc plasgraph2 37
gplascc mlplasmids 50
mobrecon rfplasmid 1
mobrecon plasclass 2
mobrecon plasgraph2 2
mobrecon mlplasmids 2
plasbinflow rfplasmid 1
plasbinflow plasclass 3
plasbinflow plasgraph2 2
plasbinflow mlplasmids 2
> grep -c SUCCESS report_plasmerge_run1_20250811.txt
7835
> grep -c ERROR report_plasmerge_run1_20250811.txt
149
```
Overall 149 experiments did not complete.  
Report in `report_plasmerge_run1_20250811.txt` and samples to re-run in `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.csv` with one line per pair `(binning method,classification method)`.
Failed samples listed in `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt`
```
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt
50 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt
```
Overall, 50 samples had at least one experiment that did not work (file `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt`) so all samples that did not work were in part due to `gplascc+mlplasmids`.

##  PlasEval (1)
`2025-08-11`: running an updated PlasEval (branch `rev_comp_output`) where the progam exits properly when the number of iterations is reached and the output is modified to include both normalized and non-normalized scores.  
Parameters: `min_len=100, alpha=0.5, max_calls=1000000`.  
Results in directory `eval/`.  
```
sbatch run_plasmerge_gt_all.sh
sbatch run_plasmerge_gp_all.sh
sbatch run_plasmerge_pbf_all.sh
sbatch run_plasmerge_mob_all.sh
```
Finished on `2028-08-12`.  
Checking results and recording success/errors.
```
> ./check_plaseval_all.sh report_plaseval_run1_20250812.txt 500
ground_truth rfplasmid unmerged 2
ground_truth rfplasmid merged 4
ground_truth plasclass unmerged 2
ground_truth plasclass merged 2
ground_truth plasgraph2 unmerged 2
ground_truth plasgraph2 merged 3
ground_truth mlplasmids unmerged 3
ground_truth mlplasmids merged 4
gplascc rfplasmid unmerged 5
gplascc rfplasmid merged 3
gplascc plasclass unmerged 39
gplascc plasclass merged 36
gplascc plasgraph2 unmerged 38
gplascc plasgraph2 merged 37
gplascc mlplasmids unmerged 56
gplascc mlplasmids merged 51
mobrecon rfplasmid unmerged 1
mobrecon rfplasmid merged 1
mobrecon plasclass unmerged 2
mobrecon plasclass merged 2
mobrecon plasgraph2 unmerged 2
mobrecon plasgraph2 merged 2
mobrecon mlplasmids unmerged 2
mobrecon mlplasmids merged 2
plasbinflow rfplasmid unmerged 43
plasbinflow rfplasmid merged 2
plasbinflow plasclass unmerged 96
plasbinflow plasclass merged 18
plasbinflow plasgraph2 unmerged 61
plasbinflow plasgraph2 merged 12
plasbinflow mlplasmids unmerged 64
plasbinflow mlplasmids merged 4
> grep -c SUCCESS report_plaseval_run1_20250812.txt
15367
> grep -c ERROR report_plaseval_run1_20250812.txt
601
```
Overall 601 experiments did not complete either due to PlasMerge results not available or PlasEval failing.  
Report in `report_plaseval_run1_20250812.txt` and samples to re-run in `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval_errors_run1.csv` with one line per pair `(binning method,classification method,merged/unmerged)`. Failed samples listed in `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval_errors_run1.samples.txt`.  
There were a lot of failed runs of PlasEval with PlasBin-flow, or even surprisingly, with unmerged ground truth bins for which it should be trivial.  
```
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval_errors_run1.samples.txt
134 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plaseval_errors_run1.samples.txt
```
Overall 134 out of 499 samples had at least one error with either PlasMerge or PlasEval.  

## Analysis
We first collect all results (PlasEval scores) into a single CSV file `analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv` with some statistics on samples with issues in in `analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.NA.txt`.

```
> source ../../env_plaseval/bin/activate
> python analysis_utils.py csv \
  plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv \
  eval \
  analysis plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv \
  499 \
> analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.NA.txt
> cat analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.NA.txt
Number of rows: 7984
Number of rows with NA: 491
unmerged.nb_bins:       149 rows with NA
unmerged.nb_ctgs:       149 rows with NA
unmerged.len_ctgs:      149 rows with NA
unmerged.Dissimilarity: 483 rows with NA
unmerged.Extra_ctgs:    483 rows with NA
unmerged.Missing_ctgs:  483 rows with NA
unmerged.Cuts:  483 rows with NA
unmerged.Joins: 483 rows with NA
unmerged.Precision:     149 rows with NA
unmerged.Recall:        149 rows with NA
unmerged.F1:    149 rows with NA
merged.nb_bins: 149 rows with NA
merged.nb_ctgs: 149 rows with NA
merged.len_ctgs:        149 rows with NA
merged.Dissimilarity:   248 rows with NA
merged.Extra_ctgs:      248 rows with NA
merged.Missing_ctgs:    248 rows with NA
merged.Cuts:    248 rows with NA
merged.Joins:   248 rows with NA
merged.Precision:       149 rows with NA
merged.Recall:  149 rows with NA
merged.F1:      149 rows with NA
```
So out of `7984` possible combinations, only `491` miss at least one result (PlasMerge or PlasEval).
They will need to be looked at, but we can proceed with enough results.

Creating scatter and difference boxplots plots `merged` versus `unmerged` for all statistics.
```
> create_scatter_plots.sh
> create_difference_plots.sh
```
All figures are in `analysis/figures`.

##  PlasEval (2)
`2025-08-13`: running PlasEval with parameters: `min_len=100, alpha=0.0, max_calls=1000000` and results in directory `eval_0/`.  
Using `alpha=0` we will then record only the number of cuts and joins instead of the scores weighted by the lengh of the contigs.
```
> tar czvf run_plaseval_all_05.tar.gz run_plaseval_*_all.sh
> sed -i 's/ALPHA=0.5/ALPHA=0.0/g' run_plaseval_gp_all.sh
> sed  's/\/eval\//\/eval_0\//g' run_plaseval_gp_all.sh
> sed -i 's/ALPHA=0.5/ALPHA=0.0/g' run_plaseval_gt_all.sh
> sed  's/\/eval\//\/eval_0\//g' run_plaseval_gt_all.sh
> sed -i 's/ALPHA=0.5/ALPHA=0.0/g' run_plaseval_mob_all.sh
> sed  's/\/eval\//\/eval_0\//g' run_plaseval_mob_all.sh
> sed -i 's/ALPHA=0.5/ALPHA=0.0/g' run_plaseval_pbf_all.sh
> sed  's/\/eval\//\/eval_0\//g' run_plaseval_pbf_all.sh
> sbatch run_plaseval_gt_all.sh
> sbatch run_plaseval_gp_all.sh
> sbatch run_plaseval_pbf_all.sh
> sbatch run_plaseval_mob_all.sh
```