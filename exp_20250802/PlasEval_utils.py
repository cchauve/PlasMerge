"""
Functions for the analysis of PlasMerge/PlasEval experimental results
"""

import os
import sys
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import argparse

""" Plasmid binning methods """
GT = "ground_truth"
GP = "gplascc"
MOB = "mobrecon"
PBF = "plasbinflow"
BINNING = [GP, GT, MOB, PBF]

""" Classfication methods """
PLSC = "plasclass"
PLSG = "plasgraph2"
MLP = "mlplasmids"
RFP = "rfplasmid"
CLASSIFICATION = [PLSC, PLSG, MLP, RFP]

""" Bins status: either original or merged with PlasMerge """
MERGED = "merged"
UNMERGED = "unmerged"
MERGING = [UNMERGED, MERGED]

""" Functions for reading input PlasEval files """

"""
PlasEval files formats
1. PlasEval bins file format:
  TSV file with fields plasmid contig  contig_len
  plasmid contig  contig_len
  MOB_AE271       Contig:105:104.622      7231
  MOB_AE271       Contig:128:69.535       1537
  MOB_AE271       Contig:38:95.9752       5297
  MOB_AE271       Contig:73:107.564       1587
  ...

2. PlasEval dissimilarity scores file format
  TSV file with 7 rows, where ows Cuts, Joins, Extra_ctgs, Mssing_ctgs, Dissimilarity
  contain 2 values: unnormalized score, normalized score
  Total_ctg_length        171818
  Total_ctg_length_alpha  990.6417018345136
  Cuts    0.0     0.0
  Joins   0.0     0.0
  Extra_ctgs      840.9486826251998       0.848892875262463
  Missing_ctgs    0.0     0.0
  Dissimilarity   840.9486826251998       0.848892875262463

3. PlasEval precision/recall/F1 file format
  TSV fle with 5 rows, 2 header rows and 3 statistics rows 
  >Overall details
  #Overall_statistic      Unwtd_statistic Wtd_statistic
  Precision       0.16666666666666666     0.033703133272368485
  Recall  1.0     1.0
  F1      0.2857142857142857      0.06520853461220594
"""

# Expeced file suffixes
FILE_SUFFIX_BINS_KEY = "bins"      # PlasEval bins
FILE_SUFFIX_SCORES_KEY = "scores"  # PlasEval dissimlarity scores
FILE_SUFFIX_STATS_KEY = "stats"    # PlasEval precision/recall/F1
FILE_SUFFIX = {
    FILE_SUFFIX_STATS_KEY: "eval.out",
    FILE_SUFFIX_SCORES_KEY: "comp.out",
    FILE_SUFFIX_BINS_KEY: "tsv"
}
FILE_SUFFIX_KEYS = list(FILE_SUFFIX.keys())
def _get_file_path(
        sample, assembler, merged, classification, binning, out_dir, file_type
):
    """ Returns the path to a PlasEval file (bins, scores or stats)
    Input:
    - sample: (str) sample name
    - assembler: (str)
    - merged: (str) in MERGING
    - classification: (str) in CLASSIFICATION
    - binning: (str) in BINNING
    - out_dir: (str) path to directory where to look for file
    - file_type: (str) in FILE_SUFFIX_KEYS
    Output:
    - (str) path to file, None if file does not exist
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

# Keys of dictionary recording statistics about bins for a sample:
BINS_DICT_NB_BINS_KEY = "nb_bins"   # Number of bins
BINS_DICT_NB_CTGS_KEY = "nb_ctgs"   # Number of contigs in bins (accounting for multiplicy)
BINS_DICT_LEN_CTGS_KEY = "len_ctgs" # Total length of contigs in bins (ibid)
BINS_DICT_KEYS = [
    BINS_DICT_NB_BINS_KEY,
    BINS_DICT_NB_CTGS_KEY,
    BINS_DICT_LEN_CTGS_KEY
]
def _read_PlasEval_bins_stats(in_bins_file):    
    """ Creates a dictionary with bins statistics from a PlasEval bins file
    Input:
    - in_bins_file: path to plasmid bins file in PlasEval bins format
    Output:
    - dict(k in BINS_DICT_KEYS: value (int))
    Exception:
    - if in_bins_file does not exist, value = np.nan
    """
    if in_bins_file is None:
        return {k: np.nan for k in BINS_DICT_KEYS}
    bins_df = pd.read_csv(in_bins_file, sep="\t", header=0)
    bins_dict = {k: 0 for k in BINS_DICT_KEYS}
    bins_list = []
    for index,row in bins_df.iterrows():
        bin_id = row["plasmid"]
        if bin_id not in bins_list:
            bins_dict[BINS_DICT_NB_BINS_KEY] += 1
            bins_list.append(bin_id)
        bins_dict[BINS_DICT_NB_CTGS_KEY] += 1
        bins_dict[BINS_DICT_LEN_CTGS_KEY] += int(row["contig_len"])
    return bins_dict

# Keys of dictionary recording PlasEval dissimilarity scores
SCORES_DICT_DISSIMILARITY_KEY = "Dissimilarity" # Dissimilarity score
SCORES_DICT_EXTRA_CTGS_KEY = "Extra_ctgs"       # Extra contigs score
SCORES_DICT_MISSING_CTGS_KEY = "Missing_ctgs"   # Missing contigs score
SCORES_DICT_CUTS_KEY = "Cuts"                   # Cuts score
SCORES_DICT_JOINS_KEY = "Joins"                 # Joins score
SCORES_DICT_KEYS = [
    SCORES_DICT_DISSIMILARITY_KEY,
    SCORES_DICT_EXTRA_CTGS_KEY,
    SCORES_DICT_MISSING_CTGS_KEY,
    SCORES_DICT_CUTS_KEY,
    SCORES_DICT_JOINS_KEY
]
def _read_PlasEval_scores(in_scores_file, normalized=True):
    """ Creates a dictionary with dissimilarity scores from a PlasEval dissimilarity file
    Input:
    - in_scores_file: path to a PlasEval dissimilarity scores file
    - normalized: boolean indicating if normalized scores are read
    Output: 
    - dict(k in SCORES_DICT_KEYS: value (float))
    Exception:
    - if in_scores_file does not exist, value = np.nan
    """
    if in_scores_file is None:
        return {k: np.nan for k in SCORES_DICT_KEYS}
    scores_dict = {}
    score_idx = {True: 2, False: 1}[normalized]
    with open(in_scores_file) as in_file:
        for line in in_file:
            _line = line.rstrip().split("\t")
            score_key = _line[0]
            if score_key in SCORES_DICT_KEYS:
                scores_dict[score_key] = float(_line[score_idx])
    return scores_dict

# Keys of dictionary recording PlasEval precision/recall/F1
STATS_DICT_PRECISION_KEY = "Precision" # Precision
STATS_DICT_RECALL_KEY = "Recall"       # Recall
STATS_DICT_F1_KEY = "F1"               # F1
STATS_DICT_KEYS = [
    STATS_DICT_PRECISION_KEY,
    STATS_DICT_RECALL_KEY,
    STATS_DICT_F1_KEY
]
def _read_PlasEval_stats(in_stats_file, weighted=True):
    """ Creates a dictionary with accuracy stats from a PlasEval stats file
    Input:
    - in_stats_file: pah to a PlasEval precision/recall/F1 file
    - weighted: boolean indicating if weighted stats are read
    Output: 
    - dict(k in STATS_DICT_KEYS: value (float))
    Exception:
    - if in_stats_file does not exist, value = np.nan
    """
    if in_stats_file is None:
        return {k: np.nan for k in STATS_DICT_KEYS}
    stats_dict = {}
    stat_idx = {True: 2, False: 1}[weighted]
    with open(in_stats_file) as in_file:
        for line in in_file:
            _line = line.rstrip().split("\t")
            stat_key = _line[0]
            if stat_key in STATS_DICT_KEYS:
                stats_dict[stat_key] = float(_line[stat_idx])
    return stats_dict

""" Functions for aggregating PlasEval files data from several samples in a CSV file """

SAMPLE_KEY = "sample"
ASSEMBLER_KEY = "assembler"
BINNING_KEY = "binning"
CLASSIFICATION_KEY = "classification"
WEIGHTED_KEY = "w"
UNWEIGHTED_KEY = "u"
NORMALIZED_KEY = "n"
UNNORMALIZED_KEY = "u"
WEIGHTHED = [UNWEIGHTED_KEY, WEIGHTED_KEY]
NORMALIZED = [UNNORMALIZED_KEY, NORMALIZED_KEY]

def _data_key(in_key, in_exp, in_wn=None):
    """ Returns the key for a statistic
    Input:
    - in_key: key in BINS_DICT_KEYS+SCORES_DICT_KEYS+STATS_DICT_KEYS
    - in_exp: in MERGING
    - in_wn: in WEIGHTHED+NORMALIZED or None (ignored if in_key in BINS_DICT_KEYS)
    Output:
    - key for the data
    """
    if in_key in BINS_DICT_KEYS or in_wn is None:
        return f"{in_exp}.{in_key}"
    else:
        return f"{in_exp}.{in_wn}.{in_key}"

def _read_sample_data(sample, assembler, classification, binning, plaseval_results_dir):
    """ Creates a dictionary with all data about a sample from PlasEval files
    Input:
    - sample: sample name
    - assembler: in [UNICYCLER, SKESA]
    - classification: in CLASSIFICATION
    - binning: in BINNING
    - plaseval_results_dir: directory where to look for all samples PlasEval results
    Output:
    - dict(
        k in [SAMPLE_KEY, ASSEMBLER_KEY, classification, binning] \
             + [_data_key(a,b,None) for a in BINS_DICT_KEYS for b in MERGING] \
             + [_data_key(a,b,c) for a in STATS_DICT_KEYS for b in MERGING for c in WEIGHTHED] \
             + [_data_key(a,b,c) for a in STATS_DICT_KEYS for b in MERGING for c in NORMALIZED]:
        value (str,str,str,str,float, ...)
      )
    """
    # Initialization of the dictionary
    sample_data = {
        SAMPLE_KEY: sample,
        ASSEMBLER_KEY: assembler,
        CLASSIFICATION_KEY: classification,
        BINNING_KEY: binning
    }
    def _read_data(in_dict, in_wn):
        for k,v in in_dict.items():
            sample_data[_data_key(k,merged,in_wn)] = v
    for merged in MERGING:
        # Reading bins data
        bins_file = _get_file_path(
            sample, assembler, merged, classification, binning, plaseval_results_dir, FILE_SUFFIX_BINS_KEY
        )
        bins_data = _read_PlasEval_bins_stats(bins_file)
        _read_data(bins_data, None)
        # Reading dissimilarity scores
        scores_file = _get_file_path(
            sample, assembler, merged, classification, binning, plaseval_results_dir, FILE_SUFFIX_SCORES_KEY
        )
        scores_data = _read_PlasEval_scores(scores_file, normalized=True)
        _read_data(scores_data, NORMALIZED_KEY)
        scores_data = _read_PlasEval_scores(scores_file, normalized=False)
        _read_data(scores_data, UNNORMALIZED_KEY)
        # Reading accuracy statistics
        stats_file = _get_file_path(
            sample, assembler, merged, classification, binning, plaseval_results_dir, FILE_SUFFIX_STATS_KEY
        )
        stats_data = _read_PlasEval_stats(stats_file, weighted=True)
        _read_data(stats_data, WEIGHTED_KEY)
        stats_data = _read_PlasEval_stats(stats_file, weighted=False)
        _read_data(stats_data, UNWEIGHTED_KEY)
    return sample_data

def aggregate_results_to_csv(
        in_samples_file, in_results_dir, out_file,
        max_idx=0, verbose=False
):
    """ Reads PlasEval data from several samples and aggregate in a single CSV file
    Input:
    - in_samples_file: path to CSV file with one row per sample with expected
      columns "species_sample" (sample name), "assembler"
    - in_results_dir: path to directory where all samples results are stored
      with PlasEval for a sample expected to be in a subdirectory "sample_assembler"
    - out_file: path to the CSV file to write
    - max_id: (int) max number of rows to read in in_samples_file, if 0 all rows are read
    - verbose: (bool) if True print statistics about missing data
    """    
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
    all_data_df.to_csv(
        out_file, sep=",", header=True, index=False
    )
    
    if verbose:
        nb_rows = all_data_df.shape[0]
        print(f"Number of rows:\t{nb_rows}")
        nb_rows_with_NA = nb_rows - all_data_df.dropna().shape[0]
        print(f"Number of rows with NA:\t{nb_rows_with_NA}")
        for col in all_data_df.columns:
            print(f"{col}:\t{all_data_df[col].isna().sum()} rows with NA")

""" Plotting functions """

def _figure_title(binning, classification, columns, nb_data_points):
        if binning is None and classification is None: title = "All samples"
        elif binning is not None and classification is None: title = binning
        elif binning is None and classification is not None: title = classification
        else: title = f"{binning}+{classification}"
        return f"{title} - {columns} (n={nb_data_points})"

def _prepare_df_to_plot(in_file, binning, classification, in_cols, norm_weight):
    def _filter_df_for_combination(in_df, in_binning, in_classification):
        if in_binning is None and in_classification is None:
            out_df = in_df
        elif in_binning is not None:
            out_df = in_df.loc[in_df[BINNING_KEY]==in_binning]
        elif in_classification is not None:
            out_df = in_df.loc[in_df[CLASSIFICATION_KEY]==in_classification]
        else:
            out_df = in_df.loc[
                (in_df[BINNING_KEY]==in_binning)
                &
                (in_df[CLASSIFICATION_KEY]==in_classification)
            ]
        return out_df
    # List of columns to consider
    unmerged_cols_list = [_data_key(col,UNMERGED,norm_weight) for col in in_cols]
    merged_cols_list = [_data_key(col,MERGED,norm_weight) for col in in_cols]
    cols_list = unmerged_cols_list+merged_cols_list
    # DataFrame, without rows with NA in considered columns and with rows to plot selected
    plot_df = _filter_df_for_combination(
        pd.read_csv(
            in_file, sep=",", header=0
        ).dropna(subset=cols_list),
        binning, classification
    )
    return plot_df,unmerged_cols_list,merged_cols_list,cols_list

def create_unmerged_vs_merged_scatter_plot(
        in_results_file, out_file,
        in_cols,
        binning=None, classification=None, norm_weight=None
):
    """ Creates a scatter plot of the sum of merged valus for cols in in_cols
    versus the sum of unmerged values and saves it in a PNG file.
    If several columns are listed in in_cols, the sum of the cols is plotted
    Input:
    - in_results_file: file containing all results
    - in_cols: list of statistics to sum and plot
    - out_file: PNG file name
    - binning, classification: combination to plot (all if both None)
    - norm_weight: version of the statistic to plot, in WEIGHTED+NORMALIZED
    """
    # Creating the dataframe to plot
    plot_df,unmerged_cols_list,merged_cols_list,cols_list = _prepare_df_to_plot(
        in_results_file, binning, classification, in_cols, norm_weight
    )
    nb_data_points = plot_df.shape[0]
    # Updating the data frame to add columns with the sums to plot
    joined_cols = "+".join(in_cols)
    unmerged_cols = _data_key(joined_cols,UNMERGED,norm_weight)
    merged_cols = _data_key(joined_cols,MERGED,norm_weight)
    if len(in_cols) > 1:
        # Creating new columns with the sum of the columns in in_cols
        plot_df.loc[:, unmerged_cols] = plot_df[unmerged_cols_list].sum(axis=1)
        plot_df.loc[:, merged_cols] = plot_df[merged_cols_list].sum(axis=1)
    # Plotting
    plot_df.plot.scatter(x=unmerged_cols, y=merged_cols, grid=True)
    xy_max = max(plot_df[unmerged_cols].max(), plot_df[merged_cols].max())
    plt.plot([0,xy_max], [0,xy_max], "k-")
    plt.xlabel(unmerged_cols)
    plt.ylabel(merged_cols)
    plt.title(
        f"{_figure_title(binning, classification, in_cols, nb_data_points)}"
    )
    plt.savefig(out_file)

def create_unmerged_vs_merged_difference_violin_plot(
        in_results_file, out_file,
        in_cols,
        binning=None, classification=None, norm_weight=None
):
    """ Creates a violin plot of the difference of merged valus for cols in in_cols
    minus the sum of unmerged values for the same columns and saves it in a PNG file
    Input:
    - in_results_file: file containing all results
    - in_cols: list of statistics to sum and plot
    - out_file: PNG file name
    - binning, classification: combination to plot (all if both None)
    - norm_weight: version of the statistic to plot, in WEIGHTED+NORMALIZED
    """
    # Creating the dataframe to plot
    plot_df,unmerged_cols_list,merged_cols_list,cols_list = _prepare_df_to_plot(
        in_results_file, binning, classification, in_cols, norm_weight
    )
    nb_data_points = plot_df.shape[0]
    # Updating the data frame to add columns with the differences to plot
    diff_cols = []
    for c in range(len(merged_cols_list)):
        unmerged_col = unmerged_cols_list[c]
        merged_col = merged_cols_list[c]
        diff_col = merged_col.replace(MERGED,"diff")
        plot_df.loc[:, diff_col] = plot_df[merged_col] - plot_df[unmerged_col]
        diff_cols.append(diff_col)
    # Recording data to plot in a list of list
    plot_data = [plot_df[col].values for col in diff_cols]        
    # Plotting
    plt.violinplot(plot_data, showmeans=True, showmedians=True)
    plt.xticks([i+1 for i in range(len(diff_cols))], diff_cols, rotation=15)
    plt.title(
        f"{_figure_title(binning, classification, in_cols, nb_data_points)}"
    )
    plt.savefig(out_file)

# Command names
CSV_CMD = "csv"
SCATTER_CMD = "scatter"
DIFF_CMD = "difference"

def main(args):
    def _cmd_figure(args_cmd, fig_type, fig_function):
        columns = args_cmd.columns.split(",")
        fig_function(
            args_cmd.data_file, args_cmd.output_file,
            columns,
            binning=args_cmd.binning,
            classification=args_cmd.classification,
            norm_weight=args_cmd.version
        )
    
    if args.command == CSV_CMD:
        aggregate_results_to_csv(
            args.samples_file,
            args.input_dir,
            args.output_file,
            max_idx=args.nb_samples,
            verbose=args.verbose
        )

    elif args.command == SCATTER_CMD:
        _cmd_figure(args, "scatter", create_unmerged_vs_merged_scatter_plot)
        
    elif args.command == DIFF_CMD:
        _cmd_figure(args, "diff", create_unmerged_vs_merged_difference_violin_plot)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="analysis_utils",
        description="PlasEval analysis tools"
    )
    subparsers = parser.add_subparsers(dest="command")

    parser_csv = subparsers.add_parser(CSV_CMD, help="Record results for seveal samples in a CSV file")
    parser_csv.add_argument("samples_file", help="Samples CSV file")
    parser_csv.add_argument("input_dir", help="Directory containing samples PlasEval results")
    parser_csv.add_argument("output_file", help="Created CVS file")
    parser_csv.add_argument("-n", "--nb_samples", type=int, default=0, help="Number of samples to read")
    parser_csv.add_argument("-v", "--verbose", action="store_true", default=False, help="Verbose mode")

    parser_scatter = subparsers.add_parser(SCATTER_CMD, help="Create merged/unmerged PNG scatter plot")
    parser_scatter.add_argument("data_file", help="Data CSV file")
    parser_scatter.add_argument("output_file", help="Output PNG file")
    parser_scatter.add_argument("columns", help="Columns to plot, separated by a comma")
    parser_scatter.add_argument("-b", "--binning", default=None, help="Binning method to consider")
    parser_scatter.add_argument("-c", "--classification", default=None, help="Classification method to consider")
    parser_scatter.add_argument("-v", "--version", default="u", help="Weighted/unweighted, normalized/unnormalized: u/n/w")

    parser_diff = subparsers.add_parser(DIFF_CMD, help="Create merged/unmerged PNG difference plot")
    parser_diff.add_argument("data_file", help="Data CSV file")
    parser_diff.add_argument("output_file", help="Output PNG file")
    parser_diff.add_argument("columns", help="Columns to plot, separated by a comma")
    parser_diff.add_argument("-b", "--binning", default=None, help="Binning method to consider")
    parser_diff.add_argument("-c", "--classification", default=None, help="Classification method to consider")
    parser_diff.add_argument("-v", "--version", default="u", help="Weighted/unweighted, normalized/unnormalized: u/n/w")

    args = parser.parse_args()

    main(args)

