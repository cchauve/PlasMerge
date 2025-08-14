"""
Functions for the analysis of experimental results
"""

import os
import sys
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import product

GT = "ground_truth"
GP = "gplascc"
MOB = "mobrecon"
PBF = "plasbinflow"
BINNING = [GP, GT, MOB, PBF]
PLSC = "plasclass"
PLSG = "plasgraph2"
MLP = "mlplasmids"
RFP = "rfplasmid"
CLASSIFICATION = [PLSC, PLSG, MLP, RFP]
GT_CLASSIFICATION = CLASSIFICATION[0]
MERGED = "merged"
UNMERGED = "unmerged"
MERGING = [UNMERGED, MERGED]

BINS_DICT_NB_BINS_KEY = "nb_bins"
BINS_DICT_NB_CTGS_KEY = "nb_ctgs"
BINS_DICT_LEN_CTGS_KEY = "len_ctgs"
BINS_DICT_KEYS = [
    BINS_DICT_NB_BINS_KEY,
    BINS_DICT_NB_CTGS_KEY,
    BINS_DICT_LEN_CTGS_KEY
]
def _get_PlasEval_bins_stats(in_bins_file, sep="\t"):    
    """
    Input:
    -plasmid bins file in PlasEval format:
    TSV file with fields plasmid contig  contig_len
    Output:
    dict{
    BINS_DICT_NB_BINS_KEY: int,
    BINS_DICT_NB_CTGS_KEY: int,
    BINS_DICT_LEN_CTGS_KEY: int
    }
    - sep: separator in file
    """
    if in_bins_file is None:
        return {k: np.nan for k in BINS_DICT_KEYS}
    bins_dict = {k: 0 for k in BINS_DICT_KEYS}
    bins_list = []
    bins_df = pd.read_csv(in_bins_file, sep=sep, header=0)
    for index,row in bins_df.iterrows():
        bin_id = row["plasmid"]
        if bin_id not in bins_list:
            bins_dict[BINS_DICT_NB_BINS_KEY] += 1
            bins_list.append(bin_id)
        bins_dict[BINS_DICT_NB_CTGS_KEY] += 1
        bins_dict[BINS_DICT_LEN_CTGS_KEY] += int(row["contig_len"])
    return bins_dict

SCORES_DICT_DISSIMILARITY_KEY = "Dissimilarity"
SCORES_DICT_EXTRA_CTGS_KEY = "Extra_ctgs"
SCORES_DICT_MISSING_CTGS_KEY = "Missing_ctgs"
SCORES_DICT_CUTS_KEY = "Cuts"
SCORES_DICT_JOINS_KEY = "Joins"
SCORES_DICT_KEYS = [
    SCORES_DICT_DISSIMILARITY_KEY,
    SCORES_DICT_EXTRA_CTGS_KEY,
    SCORES_DICT_MISSING_CTGS_KEY,
    SCORES_DICT_CUTS_KEY,
    SCORES_DICT_JOINS_KEY
]
def _get_PlasEval_scores(in_scores_file, normalized=True, sep="\t"):
    """
    Input:
    - PlasEval scores file in format
    Total_ctg_length        171818
    Total_ctg_length_alpha  990.6417018345136
    Cuts    0.0     0.0
    Joins   0.0     0.0
    Extra_ctgs      840.9486826251998       0.848892875262463
    Missing_ctgs    0.0     0.0
    Dissimilarity   840.9486826251998       0.848892875262463
    - normalized: boolean indicating if normalized scores are read
    - sep: separator in file
    Output: 
    dict{
    SCORES_DICT_DISSIMILARITY_KEY: float
    SCORES_DICT_EXTRA_KEY: float
    SCORES_DICT_MISSING_KEY: float
    SCORES_DICT_CUTS_KEY: float
    SCORES_DICT_JOINS_KEY: float
    }
    """
    if in_scores_file is None:
        return {k: np.nan for k in SCORES_DICT_KEYS}
    scores_dict = {}
    normalized_score_idx = {True: 2, False: 1}
    with open(in_scores_file) as in_file:
        for line in in_file:
            _line = line.rstrip().split(sep)
            score_key = _line[0]
            if score_key in SCORES_DICT_KEYS:
                scores_dict[score_key] = float(_line[normalized_score_idx[normalized]])
    return scores_dict

STATS_DICT_PRECISION_KEY = "Precision"
STATS_DICT_RECALL_KEY = "Recall"
STATS_DICT_F1_KEY = "F1"
STATS_DICT_KEYS = [
    STATS_DICT_PRECISION_KEY,
    STATS_DICT_RECALL_KEY,
    STATS_DICT_F1_KEY
]
def _get_PlasEval_stats(in_stats_file, weighted=True, sep="\t"):
    """
    Input:
    - PlasEval stat file in format
    >Overall details
    #Overall_statistic      Unwtd_statistic Wtd_statistic
    Precision       0.16666666666666666     0.033703133272368485
    Recall  1.0     1.0
    F1      0.2857142857142857      0.06520853461220594
    - weighted: boolean indicating if weighted stats are read
    - sep: separator in file
    Output: 
    dict{
    STATS_DICT_PRECISION_KEY: float,
    STATS_DICT_RECALL_KEY: float,
    STATS_DICT_KEYS: float
    }
    """
    if in_stats_file is None:
        return {k: np.nan for k in STATS_DICT_KEYS}
    stats_dict = {}
    weighted_stat_idx = {True: 2, False: 1}
    with open(in_stats_file) as in_file:
        for line in in_file:
            _line = line.rstrip().split(sep)
            stat_key = _line[0]
            if stat_key in STATS_DICT_KEYS:
                stats_dict[stat_key] = float(_line[weighted_stat_idx[weighted]])
    return stats_dict

FILE_SUFFIX_STATS_KEY = "stats"
FILE_SUFFIX_SCORES_KEY = "scores"
FILE_SUFFIX_BINS_KEY = "bins"
FILE_SUFFIX = {
    FILE_SUFFIX_STATS_KEY: "eval.out",
    FILE_SUFFIX_SCORES_KEY: "comp.out",
    FILE_SUFFIX_BINS_KEY: "tsv"
}
FILE_SUFFIX_KEYS = list(FILE_SUFFIX.keys())
def _get_file_path(sample, assembler, merged, classification, binning, out_dir, file_type):
    """
    Input:
    - sample: sample name
    - assembler: in [UNICYCLER, SKESA]
    - merged: in MERGING
    - classification: in CLASSIFICATION
    - binning: in BINNING
    - out_dir: directory where to look for all samples results
    - file_type: in FILE_SUFFIX_KEYS
    Output:
    path to file
    """
    file_path =  os.path.join(
        out_dir,
        f"{sample}_{assembler}",
        f"{sample}_{assembler}_{binning}_{classification}.{merged}.{FILE_SUFFIX[file_type]}"
    )
    if os.path.exists(file_path):
        return file_path
    else:
        return None

SAMPLE_KEY = "sample"
ASSEMBLER_KEY = "assembler"
BINNING_KEY = "binning"
CLASSIFICATION_KEY = "classification" 

def _read_sample_data(sample, assembler, classification, binning, out_dir):
    """
    Input:
    - sample: sample name
    - assembler: in [UNICYCLER, SKESA]
    - classification: in CLASSIFICATION
    - binning: in BINNING
    - out_dir: directory where to look for all samples results
    Output:
    dict{
    SAMPLE_KEY, ASSEMBLER_KEY,
    classification, binning,
    merged and unmerged bins stats,
    merged and unmerged binning stats,
    merged and unmerged binning scores
    }
    """
    def _data_key(in_key, in_exp):
        return f"{in_exp}.{in_key}"

    sample_data = {
        SAMPLE_KEY: sample,
        ASSEMBLER_KEY: assembler,
        CLASSIFICATION_KEY: classification,
        BINNING_KEY: binning
    }
    for merged in MERGING:
        bins_file = _get_file_path(
            sample, assembler, merged, classification, binning, plaseval_results_dir, FILE_SUFFIX_BINS_KEY
        )                
        scores_file = _get_file_path(
            sample, assembler, merged, classification, binning, plaseval_results_dir, FILE_SUFFIX_SCORES_KEY
        )
        stats_file = _get_file_path(
            sample, assembler, merged, classification, binning, plaseval_results_dir, FILE_SUFFIX_STATS_KEY
        )
        bins_data = _get_PlasEval_bins_stats(bins_file)
        for k,v in bins_data.items():
            sample_data[_data_key(k,merged)] = v
        scores_data = _get_PlasEval_scores(scores_file)
        for k,v in scores_data.items():
            sample_data[_data_key(k,merged)] = v
        stats_data = _get_PlasEval_stats(stats_file)
        for k,v in stats_data.items():
            sample_data[_data_key(k,merged)] = v                    
    return sample_data

def aggregate_results_to_csv(in_samples_file, in_results_dir, out_file, max_idx=None, verbose=False):
    samples_df = pd.read_csv(in_samples_file, sep=",", header=0)


    all_data_dict = {}
    all_data_idx = 0

    for idx,row in samples_df.iterrows():
        if max_idx is None or idx < max_idx:
            sample = row["species_sample"]
            assembler = row["assembler"]
            for binning,classification in product(BINNING, CLASSIFICATION):
                sample_data = _read_sample_data(
                    sample, assembler, classification, binning, in_results_dir
                )
                all_data_dict[all_data_idx] = sample_data
                all_data_idx +=1
            
    all_data_df = pd.DataFrame.from_dict(all_data_dict, orient="index")

    if verbose:
        nb_rows = all_data_df.shape[0]
        print(f"Number of rows:\t{nb_rows}")
        nb_rows_with_NA = nb_rows - all_data_df.dropna().shape[0]
        print(f"Number of rows with NA:\t{nb_rows_with_NA}")
        for col in all_data_df.columns:
            if "merged" in col:
                print(f"{col}:\t{all_data_df[col].isna().sum()} rows with NA")
    
    all_data_df.to_csv(
        out_file, sep=",", header=True, index=False
    )

def create_unmerged_vs_merged_scatter_plot(
        in_results_file, in_cols, out_file, out_title,
        aggregate=True, binning=None, classification=None
):
    """
    Creates a scatter plot of the sum of unmerged.in_col1 columns (unmerged values)
    versus the sum of merged.in_col (merged values) and saves it in a PNG file
    Input:
    - in_csv_file: file containing all results
    - in_cols: list of column names
    - out_file: PNG file name
    - out_title: figure title
    - aggregate: boolean
      if True: all combination (binning,classfication) in one plot
      if False: one plot per combination (binning,classfication)
    """
    in_df = pd.read_csv(in_results_file, sep=",", header=0)
    # List of columns to consider
    unmerged_cols_list = [f"unmerged.{col}" for col in in_cols]
    merged_cols_list = [f"merged.{col}" for col in in_cols]
    # Dropping rows with NA in some of these columns
    filtered_df = in_df.dropna(subset=unmerged_cols_list+merged_cols_list)
    # Updating the data frame to add columns with the sums to plot
    if len(in_cols) > 1:
        # Creating new columns with the sum of the columns in in_cols
        joined_cols = "+".join(in_cols)
        unmerged_cols = f"unmerged.{joined_cols}"
        merged_cols = f"merged.{joined_cols}"
        filtered_df.loc[:, unmerged_cols] = filtered_df[unmerged_cols_list].sum(axis=1)
        filtered_df.loc[:, merged_cols] = filtered_df[merged_cols_list].sum(axis=1)
    else:
        # No new column needs to be created
        unmerged_cols = unmerged_cols_list[0]
        merged_cols = merged_cols_list[0]

    # Plotting
    if aggregate is True:
        plot_df = filtered_df
    else:
        plot_df =  filtered_df.loc[
            (filtered_df["binning"]==binning)
            &
            (filtered_df["classification"]==classification)
        ]
    nb_data_points = plot_df.shape[0]

    plot_df.plot.scatter(x=unmerged_cols, y=merged_cols, grid=True)
    xmax = plot_df[unmerged_cols].max()
    ymax = plot_df[merged_cols].max()
    plt.plot([0,xmax],[0,ymax],"k-")
    plt.xlabel(unmerged_cols)
    plt.ylabel(merged_cols)
    plt.title(f"{out_title} (n={nb_data_points})")
    plt.savefig(out_file)

def create_unmerged_vs_merged_difference_boxplot(
        in_results_file, out_file, out_title,
        in_cols="all",
        aggregate=True, binning=None, classification=None
):
    """
    Creates a boxplot of the difference between unmeged and merged statistics
    for statistics in and saves it in a PNG file
    Input:
    - in_csv_file: file containing all results
    - out_file: PNG file name
    - out_title: figure title
    - in_cols: list of column names (default: all)
    - aggregate: boolean
      if True: all combination (binning,classfication) in one plot
      if False: one plot per combination (binning,classfication)
    """
    in_df = pd.read_csv(in_results_file, sep=",", header=0)
    # List of columns to consider
    unmerged_cols_list = [f"unmerged.{col}" for col in in_cols]
    merged_cols_list = [f"merged.{col}" for col in in_cols]
    # Dropping rows with NA in some of these columns
    filtered_df = in_df.dropna(subset=unmerged_cols_list+merged_cols_list)
    # Updating the data frame to add columns with the differences to plot
    if in_cols == "all":
        cols = BINS_DICT_KEYS+STATS_DICT_KEYS+SCORES_DICT_KEYS
    else:
        cols = in_cols
    diff_cols = [f"diff.{col}" for col in cols]
    for col in cols:
        unmerged_col = f"unmerged.{col}"
        merged_col = f"merged.{col}"
        diff_col = f"diff.{col}"
        filtered_df.loc[:, diff_col] = filtered_df[merged_col] - filtered_df[unmerged_col]

    # Plotting
    if aggregate is True:
        plot_df = filtered_df[diff_cols]
    else:
        plot_df_aux =  filtered_df.loc[
            (filtered_df["binning"]==binning)
            &
            (filtered_df["classification"]==classification)
        ]
        plot_df = plot_df_aux[diff_cols]
    nb_data_points = plot_df.shape[0]

    plot_df.plot.box(grid=True, rot=15)
    plt.xlabel(diff_cols)
    plt.title(f"{out_title} (n={nb_data_points})")
    plt.savefig(out_file)
            
if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "csv":
        samples_file = sys.argv[2]         #"plasmids_benchmarking_2025-08-02_data.filtered.randomized.csv"
        plaseval_results_dir = sys.argv[3] #"eval"
        out_dir = sys.argv[4]              #"analysis"
        out_file = sys.argv[5]             #"plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv"
        nb_samples = int(sys.argv[6])      # 499
        aggregate_results_to_csv(
            samples_file,
            plaseval_results_dir,
            os.path.join(out_dir, out_file),
            max_idx=nb_samples,
            verbose=True
        )
        
    elif cmd == "scatter_aggregated":
        csv_file = sys.argv[2]            #"analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv"
        out_dir = sys.argv[3]             #"analysis/figures"
        columns = sys.argv[4].split(",")  #"Dissimilarity"
        title1 = "+".join(columns)        
        title = f"All samples {title1}"
        out_file_name = "_".join(columns)
        out_file = os.path.join(out_dir, f"scatter_{out_file_name}_aggregated.png")
        create_unmerged_vs_merged_scatter_plot(
            csv_file, columns, out_file, title, aggregate=True
        )

    elif cmd == "scatter_combination":
        csv_file = sys.argv[2]            #"analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv"
        out_dir = sys.argv[3]             #"analysis/figures"
        columns = sys.argv[4].split(",")  #"Dissimilarity"
        binning = sys.argv[5]
        classification = sys.argv[6]
        title1 = "+".join(columns)        
        title = f"{binning}+{classification} {title1}"
        out_file_name = "_".join(columns)
        out_file = os.path.join(out_dir, f"scatter_{out_file_name}_{binning}_{classification}.png")
        create_unmerged_vs_merged_scatter_plot(
            csv_file, columns, out_file, title,
            aggregate=False, binning=binning, classification=classification
        )
        
    elif cmd == "diff_aggregated":
        csv_file = sys.argv[2]            #"analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv"
        out_dir = sys.argv[3]             #"analysis/figures"
        columns = sys.argv[4].split(",")  #"Dissimilarity,Cuts,Joins" or "all"
        title = f"Merged-unmerged - all samples"
        out_file_name = "_".join(columns)
        out_file = os.path.join(out_dir, f"diff_{out_file_name}_aggregated.png")
        create_unmerged_vs_merged_difference_boxplot(
            csv_file, out_file, title, in_cols=columns, aggregate=True
        )

    elif cmd == "diff_combination":
        csv_file = sys.argv[2]            #"analysis/plasmids_benchmarking_2025-08-02_data.filtered.randomized.results.csv"
        out_dir = sys.argv[3]             #"analysis/figures"
        columns = sys.argv[4].split(",")  #"Dissimilarity,Cuts,Joins" or "all"
        binning = sys.argv[5]
        classification = sys.argv[6]
        title = f"erged-unmerged - {binning}+{classification}"
        out_file_name = "_".join(columns)
        out_file = os.path.join(out_dir, f"diff_{out_file_name}_{binning}_{classification}.png")
        create_unmerged_vs_merged_difference_boxplot(
            csv_file, out_file, title, in_cols=columns,
            aggregate=False, binning=binning, classification=classification
        )

        
        
