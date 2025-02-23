#!/usr/bin/env python

"""
PlasMerge main script
"""

import os
import argparse

from PBF_utils import DEFAULT_SCORE_OFFSET
from io_utils import (
    UNICYCLER_TAG,
    SKESA_TAG
)
from merging import (
    merging_all_pairs,
    merge_sample
)

DEFAULT_RESIDUAL_RD_THRESHOLD = 0.05
DEFAULT_MERGE_SCORE_THRESHOLD = 0

def parse_arguments():
    description = 'PlasMerge: A tool to merge plasmid bins'

    #Parsing arguments
    parser = argparse.ArgumentParser(description=description)

    #Input
    pbm_input = parser.add_argument_group('Input')
    pbm_input.add_argument("-s", "--sample", help="Sample name")
    pbm_input.add_argument("-a", "--assembly", help="Path to assembly graph file")
    pbm_input.add_argument("-t", "--assembler", default=UNICYCLER_TAG, help=f"Assembler ({UNICYCLER_TAG}(default)/{SKESA_TAG})")    
    pbm_input.add_argument("-p", "--plasmid_scores", help="Path to plasmid score file")
    pbm_input.add_argument("-b", "--plasmid_bins", help="Path to plasmid bins file")
    pbm_input.add_argument("-r", "--source", help="Plasmid bins source")
    pbm_input.add_argument("-g", "--gc_intervals", default=None, help="GC intervals file")
    pbm_input.add_argument("-o", "--offset", type=float, default=DEFAULT_SCORE_OFFSET, help="Offset for the plasmid score term")
    pbm_input.add_argument("-rt", "--rd_threshold", type=float, default=DEFAULT_RESIDUAL_RD_THRESHOLD, help="drop contigs in MILP if residual read depth < rd_threshold")
    pbm_input.add_argument("-mt", "--merge_threshold", type=float, default=0, help="threshold on pairwise score for merging")
    pbm_input.add_argument("-l", "--scoring", action="store_true", help="flag for whether scoring is performed")
    
    #Output
    pbm_output = parser.add_argument_group('Output')
    pbm_output.add_argument("-d", "--out_dir", help="Path to Gurobi output directory")
    pbm_output.add_argument("-os", "--score_file", help="Path to output TSV file (or existing TSV if --scoring is passed)")
    pbm_output.add_argument("-om", "--merge_file", help="Path to output merger file")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    if args.scoring:
        merging_all_pairs(
            args.sample,
            args.assembly,
            args.assembler,
            args.plasmid_scores,
            args.gc_intervals,
            args.plasmid_bins,
            args.source,
            args.out_dir,
            args.score_file,
            args.rd_threshold
        )

    merge_sample(
        args.assembly,
        args.assembler,
        args.plasmid_scores,
        args.gc_intervals,
        args.plasmid_bins,
        args.source,
        args.score_file,
        args.merge_file,
        score_threshold=args.merge_threshold,
    )
