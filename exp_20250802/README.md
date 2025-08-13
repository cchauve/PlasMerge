# PlasMege experiments

## Data files
```
wget https://github.com/acme92/benchmarking/blob/benchmarking_paths/plasmerge_experiments_paths.csv
cp plasmerge_experiments_paths.csv plasmids_benchmarking_2025-08-02_data.csv
head -1 plasmids_benchmarking_2025-08-02_data.csv | sed 's/,/ /g' | awk '{for(i=1;i<=NF;i++){printf("%d. %s\n",i,$i)}}' > plasmids_benchmarking_2025-08-02_data.txt
```

## Checking data
Reverting slurm script `check_data.sh` to handle all data.
```
sbatch check_data.sh
```
Checking results
```
> wc -l plasmids_benchmarking_2025-08-02_data.csv
2483
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.csv
1234
> cut -f 2 -d " " plasmids_benchmarking_2025-08-02_data.errors.txt | grep SKESA | sort -u | wc -l
675
> cut -f 2 -d " " plasmids_benchmarking_2025-08-02_data.errors.txt | grep UNICYCLER | sort -u | wc -l
574
> cut -f 3 -d " " plasmids_benchmarking_2025-08-02_data.errors.txt | sort | awk 'BEGIN{TXT="";NB=0} {if($1!=TXT){printf("%s %d\n",TXT,NB);TXT=$1;} else{NB++;}}END{printf("%s %d\n",TXT,NB);}'
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

## Randomizing input
```
> head -1 plasmids_benchmarking_2025-08-02_data.filtered.csv > plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
> grep -v species_sample plasmids_benchmarking_2025-08-02_data.filtered.csv | shuf >> plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
1234 plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv
```

## PlasMerge
Updating the slurm script `run_convert_all.sh` to process 1234 randomized samples and convert the data in PlasMerge format.
```
sbatch run_convert_all.sh
```
No error

Updating the slurm scripts `run_plasmerge_*_all.sh` to process 500 randomized samples.
```
sbatch run_plasmerge_gt_all.sh
sbatch run_plasmerge_gp_all.sh
sbatch run_plasmerge_pbf_all.sh
sbatch run_plasmerge_mob_all.sh
```
Redone on 2025-08-09 due to inconsistency in output directories between the scripts.
Finished on 2028-08-11.
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
Report in `report_plasmerge_run1_20250811.txt` and samples to re-run in `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.csv` with one line per pair `(binning method,classification method)`. Failed samples listed in `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt`
```
> wc -l plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt
50 plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt
```
Overall, 50 samples had at least one experiment that did not work (file `plasmids_benchmarking_2025-08-02_data.filtered.randomized.plasmerge_errors_run1.samples.txt`) so all samples that did not work were in part due to `gplascc+mlplasmids`.

##  PlasEval
`2025-08-11`: running an updated PlasEval (branch `rev_comp_output`) where the progam exits properly when the number of iterations is reached and the output is modified to iclud both normalized and non-normalized scores.
Parameters: `alpha=0.5, max_calls=1000000`.
```
sbatch run_plasmerge_gt_all.sh
sbatch run_plasmerge_gp_all.sh
sbatch run_plasmerge_pbf_all.sh
sbatch run_plasmerge_mob_all.sh
```
Finished on 2028-08-12.
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
Overall 134 samples had at least one error with either PlasMerge or PlasEval.

## TODO

Collect all PlasEval results in a csv file, one line per (sample,assembler,classifiation,binning) with unmerged and merged results, full score and score components, both normalized and not normalized.
In what follows, focus on normalized scores.
Create one scatter plot per (classification,binning) with merged score versus unmerged score.
Create two charts per (merged/unmerged) with one violin plot of PlasEval score for each (classfication/binning).
Repeat above but per component of the score.
Repeat above but with on chart showing the difference merged-unmerged.
Repeat all above with not normalized scores.
