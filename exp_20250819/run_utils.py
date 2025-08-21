""" Scripts for running/checking PlasMerge and PlasEval """

import os
import sys
import argparse
import csv
import subprocess
from itertools import product

""" Plasmid binning methods """
GT = "ground_truth"
GP = "gplascc"
MOB = "mobrecon"
PBF = "plasbinflow"
BINNING = [PBF, GP, MOB, GT]

""" Classfication methods """
PLSC = "plasclass"
PLSG = "plasgraph2"
MLP = "mlplasmids"
RFP = "rfplasmid"
CLASSIFICATION = [PLSC, PLSG, MLP, RFP]

""" Data access keys """
BINNING_MAP = {PBF: "pbf", GP: "gp", GT: "gt", MOB: "mob"}
CLASSIFICATION_MAP = {PLSC: "plcl", PLSG: "plgr", MLP: "mlpl", RFP: "rfpl"}
SAMPLE = "species_sample"
ASSEMBLER = "assembler"
GFA = "gfa_100_gz"
GC_PROBS = "gc_probs"
GC_BINS="gc_bins_file"
PLASMERGE="plasmerge"
PLASEVAL="plaseval"
MERGED="merged"
UNMERGED="unmerged"
MERGING_STATUS = [MERGED,UNMERGED]
SCORES="merging_scores"
SCORES_OUT="merging_scores.tsv"
EVAL_CMD="eval"
COMP_CMD="comp"
EVAL_OUT=f"{EVAL_CMD}.out"
COMP_OUT=f"{COMP_CMD}.out"
COMP_LOG=f"{COMP_CMD}.log"

""" Directories """
OUT_DIR="out_dir"

""" Plasval default """
MIN_LEN=100
MAX_CALLS=1000000
ALPHA=0.5

""" Generic functions """

class MissingFileError(Exception):
    pass

class EmptyFileError(Exception):
    pass

def _check_files(in_files, in_msg="", exit_if_pbm=True):
    outcome = True
    for in_file in in_files:
        try:
            if not os.path.exists(in_file):
                raise MissingFileError(f"file {in_file} is missing")
            elif os.path.getsize(in_file) == 0:
                raise EmptyFileError(f"file {in_file} is empty")
            else:
                print(f"LOG: {in_msg} file {in_file} OK")
        except (MissingFileError, EmptyFileError) as e:
            print(f"ERROR: {in_msg} {e.args[0]}")
            if exit_if_pbm:
                sys.exit(1)
            else:
                outcome = False
    return outcome

# PlasMerge/PlasEval merged bins file
def _bins_file_path(sample_id, assembler, binning, classification, out_dir, merged, method):
    return os.path.join(
        out_dir,
        f"{sample_id}_{assembler}.{binning}_{classification}.{method}.{merged}.tsv"
    )
# PlasMerge merging scores file
def _merging_scores_file_path(sample_id, assembler, binning, classification, out_dir):
    return os.path.join(
        out_dir,
        f"{sample_id}_{assembler}.{binning}_{classification}.{PLASMERGE}.{SCORES_OUT}"
    )
# PlasEval files
def _plaseval_file(sample_id, assembler, binning, classification, out_dir, merged, alpha, suffix):
    if len(alpha) > 0:
        suffix2 = f"_{alpha}"
    else:
        suffix2 = ""
    return os.path.join(
        out_dir,
        f"{sample_id}_{assembler}.{binning}_{classification}.{merged}.{suffix}{suffix2}"
    )

# Reading data in samples data file, for a given sample
def _read_samples_data(in_data_file, in_sample_idx, in_alpha, out_dir):
    """
    in_alpha: (str): list of alpha values comma-separated
    """
    _check_files([in_data_file])
    all_samples_data_dict = {}
    with open(in_data_file, newline='') as csvfile:
        data_lines = list(csv.DictReader(csvfile))        
        for sample_idx in in_sample_idx:
            sample_data = data_lines[sample_idx-1]
            sample_data_dict = {}
            # Reading data in CSV samples file
            sample_data_dict[SAMPLE] = sample_data[SAMPLE]
            sample_data_dict[ASSEMBLER] = sample_data[ASSEMBLER]
            sample_data_dict[GFA] = sample_data[GFA]
            sample_data_dict[GC_PROBS] = sample_data[GC_PROBS]
            for classification in CLASSIFICATION:
                sample_data_dict[classification] = sample_data[
                    f"input_pbf_{CLASSIFICATION_MAP[classification]}"
                ]
            for (binning,classification) in product(BINNING, CLASSIFICATION):
                if binning == MOB: sample_data_key = "mob_bins"
                elif binning == GT: sample_data_key = "gt_bins"
                else: sample_data_key = f"{BINNING_MAP[binning]}_{CLASSIFICATION_MAP[classification]}_bins"
                sample_data_dict[(binning,classification)] = sample_data[sample_data_key]
            # Output directory
            sample_data_dict[OUT_DIR] = os.path.join(
                out_dir,
                f"{sample_data_dict[SAMPLE]}_{sample_data_dict[ASSEMBLER]}"
            )
            # PlasMerge files in PlasMerge format
            for (binning,classification) in product(BINNING, CLASSIFICATION):
                # Merging scores file
                key = (PLASMERGE,binning,classification,SCORES)
                sample_data_dict[key] =  _merging_scores_file_path(
                    sample_data_dict[SAMPLE], sample_data_dict[ASSEMBLER],
                    binning, classification, sample_data_dict[OUT_DIR]
                )
                # Bins files
                for merged in MERGING_STATUS:
                    key = (PLASMERGE,binning,classification,merged)
                    sample_data_dict[key] = _bins_file_path(
                        sample_data_dict[SAMPLE], sample_data_dict[ASSEMBLER],
                        binning, classification, sample_data_dict[OUT_DIR],
                        merged, PLASMERGE
                    )
            # PlasEval files in PlasEval format
            for (binning,classification,merged) in product(BINNING,CLASSIFICATION,MERGING_STATUS):
                key = (PLASEVAL,binning,classification,merged)
                # Bins files
                sample_data_dict[key] = _bins_file_path(
                sample_data_dict[SAMPLE], sample_data_dict[ASSEMBLER],
                    binning, classification, sample_data_dict[OUT_DIR],
                    merged, PLASEVAL
                )
                # Eval files
                key = (PLASEVAL,binning,classification,merged,EVAL_OUT)
                sample_data_dict[key] = _plaseval_file(
                    sample_data_dict[SAMPLE], sample_data_dict[ASSEMBLER],
                    binning, classification, sample_data_dict[OUT_DIR],
                    merged, "", EVAL_OUT
                )
                if len(in_alpha) > 0:
                    # Comp files
                    for (alpha,comp_suffix) in product(in_alpha.split(","),[COMP_OUT,COMP_LOG]):
                        key = (PLASEVAL,binning,classification,merged,comp_suffix,float(alpha))
                        sample_data_dict[key] = _plaseval_file(
                            sample_data_dict[SAMPLE], sample_data_dict[ASSEMBLER],
                            binning, classification, sample_data_dict[OUT_DIR],
                            merged, alpha, comp_suffix
                        )
            all_samples_data_dict[sample_idx] = sample_data_dict
    return all_samples_data_dict

def _prefix(sample_data_dict, binning, classification):
    return f"{sample_data_dict[SAMPLE]}_{sample_data_dict[ASSEMBLER]}.{binning}_{classification}"

""" PlasMerge functions """

def _create_convert_bins_to_plasmerge_command(sample_data_dict, args):
    convert_cmd = [
        f"python{args.python_version}",
        os.path.join(args.plasmerge_path, "convert_utils.py"),
        sample_data_dict[(args.binning,args.classification)],
        sample_data_dict[(PLASMERGE,args.binning,args.classification,UNMERGED)],
        "gt"
    ]
    in_file = sample_data_dict[(args.binning,args.classification)]
    out_file = sample_data_dict[(PLASMERGE,args.binning,args.classification,UNMERGED)]
    return convert_cmd,in_file,out_file
    
def convert_bins_to_plasmerge(sample_data_dict, args):
    """ Convert a bins file into PlasMerge format
    """
    convert_cmd,in_file,out_file = _create_convert_bins_to_plasmerge_command(sample_data_dict, args)
    # Checking input files
    if os.path.exists(out_file):
        print(f"WARNING: Converting unmerged bins - file {out_file} does exist already")
        return
    _check_files([in_file], in_msg="Converting unmerged bins -")
    try:        
        convert_cmd_str = " ".join(convert_cmd)
        print(f"LOG: Converting bins - [{convert_cmd_str}]")
        result = subprocess.run(convert_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Converting bins - [{convert_cmd_str}] failed with return code {e.returncode}")
        sys.stderr.write(e.stderr)
        sys.exit(1)
    _check_files([out_file], in_msg="Converting unmerged bins -")

def _create_plasmerge_command(sample_data_dict, args):
    gurobi_out_dir = os.path.join(
        sample_data_dict[OUT_DIR],
        _prefix(sample_data_dict, args.binning, args.classification)
    )
    os.makedirs(gurobi_out_dir, exist_ok=True)
    # Command
    plasmerge_cmd = [
        f"python{args.python_version}",
        os.path.join(args.plasmerge_path, "plasmerge.py"),
        "-s", sample_data_dict[SAMPLE],
        "-t", sample_data_dict[ASSEMBLER].lower(),
        "-r", BINNING_MAP[args.binning],
        "-a", sample_data_dict[GFA],
        "-p", sample_data_dict[args.classification],
        "-b", sample_data_dict[(PLASMERGE,args.binning,args.classification,UNMERGED)],
        "-l",
        "-d", gurobi_out_dir,
        "-os", sample_data_dict[(PLASMERGE,args.binning,args.classification,SCORES)],
        "-om", sample_data_dict[(PLASMERGE,args.binning,args.classification,MERGED)],
        "-g", args.gc_bins_file
    ]
    in_files = [
        sample_data_dict[GFA],
        sample_data_dict[args.classification],
        sample_data_dict[(PLASMERGE,args.binning,args.classification,UNMERGED)],
        args.gc_bins_file
    ]
    out_files = [
        sample_data_dict[(PLASMERGE,args.binning,args.classification,MERGED)],
        sample_data_dict[(PLASMERGE,args.binning,args.classification,SCORES)]
    ]
    return plasmerge_cmd,in_files,out_files

def run_plasmerge(args):
    """ Run PlasMerge on a single sample for a combination (binning,classification)
    """
    # Reading data for the sample to process
    sample_data_dict = _read_samples_data(
        args.data_file, [args.sample_idx], "", args.output_dir
    )[args.sample_idx]
    # Converting unmerged bins into PlasMerge format
    os.makedirs(sample_data_dict[OUT_DIR], exist_ok=True)
    convert_bins_to_plasmerge(sample_data_dict, args)    
    # Running PlasMerge
    plasmerge_cmd,in_files,out_files = _create_plasmerge_command(sample_data_dict, args)
    _check_files(in_file, in_msg="Running PlasMerge -")
    try:
        plasmerge_cmd_str = " ".join(plasmerge_cmd)
        print(f"LOG: Running PlasMerge - [{plasmerge_cmd_str}]")
        result = subprocess.run(plasmerge_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Running PlasMerge - [{plasmerge_cmd_str}] failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    _check_files(out_file, in_msg="Running PlasMerge -")

def check_plasmerge(args):
    """ Check PlasMerge input/output on a set of samples for all combinations (binning,classification)
    """
    samples_id_range = range(1,args.nb_samples+1)
    all_samples_data_dict = _read_samples_data(
        args.data_file, samples_id_range, "", args.output_dir
    )
    id_with_problems = []
    for (sample_idx,binning,classification) in product(
            samples_id_range,args.binning.split(","),args.classification.split(",")
    ):
        sample_data_dict = all_samples_data_dict[sample_idx]
        files_to_check = [
            sample_data_dict[GFA],
            sample_data_dict[classification],
            sample_data_dict[(PLASMERGE,binning,classification,UNMERGED)],
            sample_data_dict[(PLASMERGE,binning,classification,MERGED)],
            sample_data_dict[(PLASMERGE,binning,classification,SCORES)]
        ]
        correct = _check_files(
            files_to_check,
            in_msg=_prefix(sample_data_dict, binning, classification),
            exit_if_pbm=False
        )
        if not correct:
            id_with_problems.append((sample_idx,binning,classification))
    with open(args.data_file) as in_file, open(args.output_file,"w") as out_file:
        in_lines = in_file.readlines()
        out_file.write(f"{in_lines[0].rstrip()},binning,classification")
        for (sample_idx,binning,classification) in id_with_problems:
            out_file.write(f"\n{in_lines[sample_idx].rstrip()},{binning},{classification}")
            
""" PlasEval functions """

def _create_convert_bins_to_plaseval_command(sample_data_dict, args):
    convert_cmd = [
        f"python{args.python_version}",
        os.path.join(args.plasmerge_path, "convert_utils.py"),
        sample_data_dict[(PLASMERGE,args.binning,args.classification,args.merged_status)],
        sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status)],
        "to_plaseval",
        sample_data_dict[GFA]
    ]
    in_files = [
        sample_data_dict[GFA],
        sample_data_dict[(PLASMERGE,args.binning,args.classification,args.merged_status)],
    ]
    out_file = sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status)]
    return convert_cmd,in_files,out_file

def convert_bins_to_plaseval(sample_data_dict, args):
    """ Convert a bins file into PlasEval format
    """
    convert_cmd,in_files,out_file = _create_convert_bins_to_plaseval_command(
        sample_data_dict, args
    )
    if os.path.exists(out_file):
        print(f"WARNING: Converting bins - file {out_file} does exist already")
        return
    _check_files(in_files, in_msg="Converting bins -")
    try:       
        convert_cmd_str = " ".join(convert_cmd)
        print(f"LOG: Converting bins - [{convert_cmd_str}]")
        result = subprocess.run(convert_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Converting bins - [{convert_cmd_str}] failed with return code {e.returncode}")
        sys.stderr.write(e.stderr)
        sys.exit(1)
    _check_files([out_file], in_msg="Converting bins -")

def _create_plaseval_eval_command(sample_data_dict, args):
    plaseval_cmd = [
        f"python{args.python_version}",
        os.path.join(args.plaseval_path, "plaseval.py"),
        EVAL_CMD,
        "--pred", sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status)],
        "--gt", sample_data_dict[(args.binning,args.classification)],
        "--out_file", sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status,EVAL_OUT)],
        "--min_len", str(args.min_len)
    ]
    in_files = [
        sample_data_dict[(args.binning,args.classification)],
        sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status)]
    ]
    out_file = sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status,EVAL_OUT)]
    return plaseval_cmd,in_files,out_file

def _run_plaseval_eval(args):
    """ Run PlasEval eval on a single sample for a combination (binning,classification)
    """
    # Reading data for the sample to process
    sample_data_dict = _read_samples_data(
        args.data_file, [args.sample_idx], "", args.output_dir
    )[args.sample_idx]
    # Converting unmerged bins into PlasMerge format
    convert_bins_to_plaseval(sample_data_dict, args)    
    # Running PlasMerge
    plaseval_cmd,in_files,out_file = _create_plaseval_eval_command(sample_data_dict, args)
    _check_files(in_files, in_msg="Running PlasEval -")
    try:
        plaseval_cmd_str = " ".join(plaseval_cmd)
        print(f"LOG: Running PlasEval - [{plaseval_cmd_str}]")
        result = subprocess.run(plaseval_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Running PlasEval - [{plaseval_cmd_str}] failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    _check_files([out_file], in_msg="Running PlasEval -")

def _create_plaseval_comp_command(sample_data_dict, args):
    plaseval_cmd = [
        f"python{args.python_version}",
        os.path.join(args.plaseval_path, "plaseval.py"),
        COMP_CMD,
        "--l", sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status)],
        "--r", sample_data_dict[(args.binning,args.classification)],
        "--out_file", sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status,COMP_OUT,args.alpha)],
        "--log_file", sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status,COMP_LOG,args.alpha)],
        "--p", str(args.alpha),
        "--min_len", str(args.min_len),
        "--max_calls", str(args.max_calls)
    ]
    in_files = [
        sample_data_dict[(args.binning,args.classification)],
        sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status)]
    ]
    out_files = [
        sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status,COMP_OUT,args.alpha)],
        sample_data_dict[(PLASEVAL,args.binning,args.classification,args.merged_status,COMP_LOG,args.alpha)]
    ]
    return plaseval_cmd,in_files,out_files

def _run_plaseval_comp(args):
    """ Run PlasEval comp on a single sample for a combination (binning,classification)
    """
    # Reading data for the sample to process
    sample_data_dict = _read_samples_data(
        args.data_file, [args.sample_idx], str(args.alpha), args.output_dir
    )[args.sample_idx]
    # Converting unmerged bins into PlasMerge format
    convert_bins_to_plaseval(sample_data_dict, args)    
    # Running PlasMerge
    plaseval_cmd,in_files,out_files = _create_plaseval_comp_command(sample_data_dict, args)
    _check_files(in_files, in_msg="Running PlasEval -")
    try:
        plaseval_cmd_str = " ".join(plaseval_cmd)
        print(f"LOG: Running PlasEval - [{plaseval_cmd_str}]")
        result = subprocess.run(plaseval_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Running PlasEval - [{plaseval_cmd_str}] failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    _check_files(out_files, in_msg="Running PlasEval -")

def run_plaseval(args):
    if args.plaseval_mode == EVAL_CMD:
        _run_plaseval_eval(args)
    elif args.plaseval_mode == COMP_CMD:
        _run_plaseval_comp(args)
    
def check_plaseval(args):
    """ Check PlasEval input/output on a set of samples for all combinations (binning,classification)
    """
    samples_id_range = range(1,args.nb_samples+1)
    all_samples_data_dict = _read_samples_data(
        args.data_file, samples_id_range, args.alpha, args.output_dir
    )
    id_with_problems = []
    for (sample_idx,binning,classification,merged_status) in product(
            samples_id_range,args.binning.split(","),args.classification.split(","),args.merged_status.split(",")
    ):
        sample_data_dict = all_samples_data_dict[sample_idx]
        # Checking PlasEval input files
        files_to_check = [
            sample_data_dict[GFA],
            sample_data_dict[(binning,classification)],
            sample_data_dict[(PLASMERGE,binning,classification,merged_status)],
            sample_data_dict[(PLASEVAL,binning,classification,merged_status)]
        ]
        input_correct = _check_files(
            files_to_check,
            in_msg=_prefix(sample_data_dict, binning, classification),
            exit_if_pbm=False
        )
        if not input_correct:
            id_with_problems.append((sample_idx,args.binning,args.classification,merged_status,"input"))        
        # Checking eval files
        files_to_check = [
            sample_data_dict[(PLASEVAL,binning,classification,merged_status,EVAL_OUT)]
        ]
        eval_correct = _check_files(
            files_to_check,
            in_msg=_prefix(sample_data_dict, binning, classification),
            exit_if_pbm=False
        )        
        if not eval_correct:
            id_with_problems.append((sample_idx,binning,classification,merged_status,EVAL_CMD))
        # Checking comp files
        for alpha in args.alpha.split(","):
            files_to_check = [
                sample_data_dict[(PLASEVAL,binning,classification,merged_status,COMP_OUT,float(alpha))],
                sample_data_dict[(PLASEVAL,binning,classification,merged_status,COMP_LOG,float(alpha))]
            ]
            comp_correct = _check_files(
                files_to_check,
                in_msg=_prefix(sample_data_dict, binning, classification),
                exit_if_pbm=False
            )
            if not comp_correct:
                id_with_problems.append((sample_idx,binning,classification,merged_status,alpha))
    with open(args.data_file) as in_file, open(args.output_file,"w") as out_file:
        in_lines = in_file.readlines()
        out_file.write(f"{in_lines[0]},binning,classification,merged,alpha")
        for (sample_idx,binning,classification,merged,alpha) in id_with_problems:
            out_file.write(f"\n{in_lines[sample_idx]},{binning},{classification},{merged},{alpha}")

""" Main function """

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="run_utils",
        description="Running PlasMege and PlasEval"
    )
    subparsers = parser.add_subparsers(dest="command")

    parser_run_plasmerge = subparsers.add_parser(PLASMERGE, help="Running PlasMerge")
    parser_run_plasmerge.add_argument("-d", "--data_file", help="Samples CSV file")
    parser_run_plasmerge.add_argument("-i", "--sample_idx", type=int, help="Index of sample to process")
    parser_run_plasmerge.add_argument("-b", "--binning", default=None, help="Binning method to consider")
    parser_run_plasmerge.add_argument("-c", "--classification", default=None, help="Classification method to consider")
    parser_run_plasmerge.add_argument("-p", "--plasmerge_path", help="Path to PlasMerge scripts")
    parser_run_plasmerge.add_argument("-gc", "--gc_bins_file", default="gc.txt", help="Path to GC bins file")    
    parser_run_plasmerge.add_argument("-o", "--output_dir", help="Directory where all results are written")
    parser_run_plasmerge.add_argument("-v", "--python_version", default="3.9", help="Python vesion")

    CHECK_PLASMERGE = f"check_{PLASMERGE}"
    parser_check_plasmerge = subparsers.add_parser(CHECK_PLASMERGE, help="Checking PlasMerge")
    parser_check_plasmerge.add_argument("-d", "--data_file", help="Samples CSV file")
    parser_check_plasmerge.add_argument("-n", "--nb_samples", type=int, help="Number of sample to check")
    parser_check_plasmerge.add_argument("-b", "--binning", default=None, help="Binning methods to consider, comma-separated")
    parser_check_plasmerge.add_argument("-c", "--classification", default=None, help="Classifications method to consider, comma-separated")
    parser_check_plasmerge.add_argument("-od", "--output_dir", help="Directory where all results are written")
    parser_check_plasmerge.add_argument("-of", "--output_file", help="File where to write samples with errors")

    parser_run_plaseval = subparsers.add_parser(PLASEVAL, help="Running PlasEval")
    parser_run_plaseval.add_argument("-d", "--data_file", help="Samples CSV file")
    parser_run_plaseval.add_argument("-i", "--sample_idx", type=int, help="Index of sample to process")
    parser_run_plaseval.add_argument("-b", "--binning", default=None, help="Binning method to consider")
    parser_run_plaseval.add_argument("-c", "--classification", default=None, help="Classification method to consider")
    parser_run_plaseval.add_argument("-pm", "--plasmerge_path", help="Path to PlasMerge scripts")
    parser_run_plaseval.add_argument("-pe", "--plaseval_path", help="Path to PlasEval scripts")
    parser_run_plaseval.add_argument("-m", "--plaseval_mode", default="eval", help="PlasEval mode (eval/comp)")
    parser_run_plaseval.add_argument("-bm", "--merged_status", default=UNMERGED, help="Merged/unmerged bins")
    parser_run_plaseval.add_argument("-a", "--alpha", default=ALPHA, type=float, help="PlasEval alpha parameter")
    parser_run_plaseval.add_argument("-ml", "--min_len", default=MIN_LEN, type=int, help="PlasEval min_len parameter")
    parser_run_plaseval.add_argument("-mc", "--max_calls", default=MAX_CALLS, type=int, help="PlasEval max_calls parameter")    
    parser_run_plaseval.add_argument("-o", "--output_dir", help="Directory where all results are written")
    parser_run_plaseval.add_argument("-v", "--python_version", default="3.9", help="Python vesion")    

    CHECK_PLASEVAL = f"check_{PLASEVAL}"
    parser_check_plaseval = subparsers.add_parser(CHECK_PLASEVAL, help="Checking PlasMerge")
    parser_check_plaseval.add_argument("-d", "--data_file", help="Samples CSV file")
    parser_check_plaseval.add_argument("-n", "--nb_samples", type=int, help="Number of sample to check")
    parser_check_plaseval.add_argument("-b", "--binning", default=None, help="Binning methods to consider, comma-separated")
    parser_check_plaseval.add_argument("-c", "--classification", default=None, help="Classifications method to consider, comma-separated")
    parser_check_plaseval.add_argument("-bm", "--merged_status", default=UNMERGED, help="Merged and/or unmerged bins, comma-separated")
    parser_check_plaseval.add_argument("-a", "--alpha", default="", help="Alpha values, comma-separated")
    parser_check_plaseval.add_argument("-ml", "--min_len", default=MIN_LEN, type=int, help="PlasEval min_len parameter")
    parser_check_plaseval.add_argument("-mc", "--max_calls", default=MAX_CALLS, type=int, help="PlasEval max_calls parameter")        
    parser_check_plaseval.add_argument("-od", "--output_dir", help="Directory where all results are written")
    parser_check_plaseval.add_argument("-of", "--output_file", help="File where to write samples with errors")

    args = parser.parse_args()

    if args.command == PLASMERGE:
        run_plasmerge(args)

    elif args.command == CHECK_PLASMERGE:
        check_plasmerge(args)

    elif args.command == PLASEVAL:
        run_plaseval(args)

    elif args.command == CHECK_PLASEVAL:
        check_plaseval(args)
