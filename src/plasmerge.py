#!/usr/bin/env python

"""
PlasMerge main script
"""

import os
import argparse

from PBF_utils import DEFAULT_SCORE_OFFSET
from merging import merging_all_pairs

DEFAULT_THRESHOLD = 0.05

def parse_arguments():
    description = 'PlasMerge: A tool to merge plasmid bins'

    #Parsing arguments
    parser = argparse.ArgumentParser(description=description)

    #Input
    pbm_input = parser.add_argument_group('Input')
    pbm_input.add_argument("-s", "--sample", help="Sample name")
    pbm_input.add_argument("-a", "--assembly", help="Path to assembly graph file")
    pbm_input.add_argument("-p", "--plasmid_scores", help="Path to plasmid score file")
    pbm_input.add_argument("-b", "--plasmid_bins", help="Path to plasmid bins file")
    pbm_input.add_argument("-r", "--source", help="Plasmid bins source")
    pbm_input.add_argument("-g", "--gc_intervals", default=None, help="GC intervals file")
    pbm_input.add_argument("-o", "--offset", type=float, default=DEFAULT_SCORE_OFFSET, help="Offset for the plasmid score term")
    pbm_input.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD, help="????")
    
    #Output
    pbm_output = parser.add_argument_group('Output')
    pbm_output.add_argument("-d", "--out_dir", help="Path to Gurobi output directory")
    pbm_output.add_argument("-f", "--out_file", help="Path to output TSV file")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    merging_all_pairs(
        args.sample,
        args.assembly,
        args.plasmid_scores,
        args.gc_intervals,
        args.plasmid_bins,
        args.source,
        args.out_dir,
        args.out_file,
        args.threshold
    )
