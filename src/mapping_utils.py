#!/usr/bin/env python

'''
Manipulating mapping contigs/genomes to define assembly ground truth
'''

import csv
import gzip
import os
import sys
import pandas as pd
import argparse
import subprocess
from collections import defaultdict

from assembly_utils import read_FASTA_len
from df_utils import (
    add_column_default,
    add_column_function,
    modify_column_function
)

# Separators for CSV/TSV files
PLASGRAPH_SEP1 = ','
PLASBIN_SEP1 = '\t'
INTERVALS_SEP = ';'
INTERVAL_SEP = ':'
# Molecule types
MOL_CHR = 'chromosome'
MOL_PLS = 'plasmid'
MOL_UNKN = 'unknown'

# Mapping contigs against a molecule (chromosome(s) or plasmid(s))

def __run_cmd_stdout(in_cmd, in_out_file, max_attempts):
    '''
    Running in_cmd (list(str)) with output in in_out_file
    Tries at most max_attempts times in case of exception raised
    '''
    with open(in_out_file, 'w') as out_file:
        for attempt in range(max_attempts):
            try:
                subprocess.run(in_cmd, check=True, stdout=out_file)
                success = True
                break
            except subprocess.CalledProcessError:
                success = False
    if not success:
        print(f'ERROR\tcreating {in_out_file}', file=sys.stderr)
    return success

def run_ctg_mol_blast(
        in_ctg_file, in_mol_file, out_tsv_file,
        in_blast_exec, in_blast_options, 
        max_attempts=5
):
    '''
    BLAST contigs to molecule.
    Input:
    - in_ctg_file: path to a FASTA file containing contigs
    - in_mol_file: path to a FASTA file containing the molecule
    - out_tsv_file: path to the TSV file that will record the results
    - in_blast_exec: path to blastn executable file
    - in_blast_options: list(str) blast options 
    - max_attempts: (int) maximum number of times BLAST is tried
    '''
    BLAST_OUT6_FORMAT = '6 qseqid qlen sseqid slen qstart qend sstart send nident sstrand length'
    blast_input = ['-query', in_ctg_file, '-subject', in_mol_file]
    blast_options = in_blast_options
    blast_out = ['-outfmt', BLAST_OUT6_FORMAT]
    blast_cmd = [in_blast_exec]+blast_input+blast_options+blast_out
    return __run_cmd_stdout(blast_cmd, out_tsv_file, max_attempts)

def run_ctg_mol_minimap2(
        in_ctg_file, in_mol_file, out_paf_file,
        in_minimap2_exec, in_minimap2_options, 
        max_attempts=5
):
    '''
    Map contigs to molecule using minimap2.
    Input:
    - in_ctg_file: path to a FASTA file containing contigs
    - in_mol_file: path to a FASTA file containing the molecule
    - out_file: path to the PA/TSVF file that will record the results
    - in_minimap2_exec: path to minimap2 executable file
    - in_minimap2_options: list(str) minimap2 options 
    - max_attempts: (int) maximum number of times BLAST is tried
    '''
    minimap2_input = [in_mol_file, in_ctg_file]
    minimap2_options = in_minimap2_options
    minimap2_cmd = [in_minimap2_exec]+minimap2_options+minimap2_input
    return __run_cmd_stdout(minimap2_cmd, out_paf_file, max_attempts)

def run_ctg_mol_mapper(
        in_mapper, 
        in_ctg_file, in_mol_file, out_file,
        in_mapper_exec, in_mapper_options=[]
):
    '''
    Map contigs to molecule using minimap2 or blast
    Input:
    - in_mapper: blast or minimap2
    - in_ctg_file: path to a FASTA file containing contigs
    - in_mol_file: path to a FASTA file containing the molecule
    - out_file: path to the file that will record the results
    - in_mapper_exec: path to mapper executable file
    - in_mapper_options: list(str) mapper options 
    '''
    subprocess.run(['gunzip', in_ctg_file], check=True)
    ctg_file = os.path.splitext(in_ctg_file)[0]
    subprocess.run(['gunzip', in_mol_file], check=True)
    mol_file = os.path.splitext(in_mol_file)[0]
    if in_mapper == 'blast':
        success = run_ctg_mol_blast(
            ctg_file, mol_file, out_file,
            in_mapper_exec, in_mapper_options
        )
    elif in_mapper == 'minimap2':
        success = run_ctg_mol_minimap2(
            ctg_file, mol_file, out_file,
            in_mapper_exec, in_mapper_options
        )    
    else:
        print(f'ERROR\t{in_mapper} must be "blast" or "minimap2"')
        sys.exit(1)
    subprocess.run(['gzip', ctg_file], check=True)
    subprocess.run(['gzip', mol_file], check=True)
    if not success:
        print(f'ERROR\tmapping {in_mapper}', file=sys.stderr)
        sys.exit(1)
        
# Reading mapping files

'''
Reading contigs to molecules (chromosome or plasmid) mapping files
The result of reading a mapping file is a pandas dataframe with the 
following columns, where each entry represents an alignment/mapping:
- ctg_id:   contig
- ctg_len:  length of contig
- mol_id:   (GenBank/RefSeq) id of the target molecule (chromosome/plasmid)
- mol_len:  length of molecule
- mol_type: in {MOL_CHR,MOL_PLS,MOL_UNKN}
- ident:    percentage of identity of the alignment defined 
  - BLAST:   nident / length
  - minimap: nmatch / alen
- ctg_[start,end]: coordinates on the contig of the alignment
- mol_[start,end]: coordinates on the molecule of the alignment
- strand:   aligment strand on molecule, in {+,-}
- length:   alignment length
  - BLAST:   length field
  - minimap: alen field

Coordinates are 1-based.
'''

def __order_coordinates(in_df, start_col, end_col):
    '''
    Ensures that the colums start_col,end_col are increasing
    '''
    for idx,row in in_df.iterrows():
        start,end = row[start_col],row[end_col]
        if start>end:
            start,end=end,start
        in_df.at[idx,start_col] = start
        in_df.at[idx,end_col] = end

# Reading a BLAST file
# Comments:
# BLAST 6 format: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# Specific required format: See BLAST_OUT6_FORMAT in run_ctg_mol_blast
# Column nident is used for the number of matches
BLAST6_COL_NAMES = [
    'ctg_id', 'ctg_len', 'mol_id', 'mol_len', 
    'ctg_start', 'ctg_end', 'mol_start', 'mol_end',
    'nident', 'strand', 'length'
]
BLAST6_COL_TYPES = {
    'ctg_id': str, 'ctg_len': int, 'mol_id': str, 'mol_len': int,
    'ctg_start': int, 'ctg_end': int, 'mol_start': int, 'mol_end': int,
    'nident': int, 'strand': str, 'length': int
}

def __read_blast(in_blast_file, mol_type=MOL_UNKN):
    '''
    Read BLAST 6 output format expected to be generated with columns
    "6 qseqid qlen sseqid slen qstart qend sstart send nident sstrand length"
    in this order
    '''
    mappings_df = pd.read_csv(
        in_blast_file, sep='\t',
        names=BLAST6_COL_NAMES,
        dtype=BLAST6_COL_TYPES
    )
    add_column_default(mappings_df, 'mol_type', default_val=mol_type)
    modify_column_function(
        mappings_df, 'strand',
        lambda x: {'plus': '+', 'minus': '-'}[x]
    )        
    add_column_function(
        mappings_df, 'ident', ['nident','length'],
        lambda x: float(x[0]/x[1])
    )
    __order_coordinates(mappings_df, 'ctg_start', 'ctg_end')
    __order_coordinates(mappings_df, 'mol_start', 'mol_end')
    return mappings_df

# Comments:
# PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md
# Columns NM to cg are SAM format tags, not all of them are always present
# Column nmatch is used for the number of matches
MINIMAP_PAF_COL_NAMES = [
    'ctg_id', 'ctg_len', 'ctg_start', 'ctg_end', 'strand', 'mol_id', 'mol_len',
    'mol_start', 'mol_end', 'nmatch', 'length'
]
MINIMAP_PAF_COL_TYPES = {
    'ctg_id': str, 'mol_id': str, 'ctg_start': int, 'ctg_end': int,
    'mol_start': int, 'mol_end': int, 'ctg_len': int, 'mol_len': int,
    'strand': str, 'nmatch': int, 'length': int
}
    
def __read_minimap2_paf(in_paf_file, mol_type=MOL_UNKN):
    mappings_df = pd.read_csv(
        in_paf_file, sep='\t',
        names=MINIMAP_PAF_COL_NAMES,
        usecols=MINIMAP_PAF_COL_NAMES,
        dtype=MINIMAP_PAF_COL_TYPES
    )
    add_column_default(mappings_df, 'mol_type', default_val=mol_type)
    add_column_function(
        mappings_df, 'ident', ['nmatch','length'],
        lambda x: float(x[0]/x[1])
    )
    # Updating coordinates to make them 1-based
    modify_column_function(mappings_df, 'ctg_start', lambda x: x+1)
    modify_column_function(mappings_df, 'mol_start', lambda x: x+1)
    __order_coordinates(mappings_df, 'ctg_start', 'ctg_end')
    __order_coordinates(mappings_df, 'mol_start', 'mol_end')
    return mappings_df   

def read_mappings_file(in_file, in_format, mol_type):
    '''
    Read a mapping file
    input:
    - in_file:   input mapping file
    - in_format: in {blast,minimap2}
    - mol_type:  in {plasmid,chromosome,unknown}
    Output:
    - dataframe as described above
    '''
    if os.path.isfile(in_file):
        if in_format == 'blast':
            return __read_blast(in_file, mol_type=mol_type)
        elif in_format == 'minimap2':
            return __read_minimap2_paf(in_file, mol_type=mol_type)
    else:
        return pd.DataFrame(
            data=None,
            columns=MINIMAP_PAF_COL_NAMES
        )
        
# Filtering mappings

def filter_low_quality_mappings(in_mappings_df, min_len, min_identity):
    '''
    Filtering mappings s.t.
    - either length on contig or on molecule < min_len 
    - or ident < min_identity
    '''
    idx_to_drop = []
    for idx,row in in_mappings_df.iterrows():        
        if row['ident'] < min_identity:
            idx_to_drop.append(idx)
        else:
            ctg_mapping_len = row['ctg_end']-row['ctg_start']+1
            mol_mapping_len = row['mol_end']-row['mol_start']+1
            if min(ctg_mapping_len,mol_mapping_len) < min_len:
                idx_to_drop.append(idx)
    in_mappings_df.drop(idx_to_drop, axis=0, inplace=True)
    
# Computing covered intervals for pairs (contig,molecule)
# A covered interval is defined by a st of mappings that cover a maximal contiguous
# segment of the considered sequence (contig or molecule) through overlaps or being
# consecutive with no gap.

def __compute_ctg_mol_intervals(ctg_mol_mappings_df, col_start, col_end):
    '''
    Input:
    - ctg_mol_mappings_df: dataframe of mappings restricted to a single pair (contig,molecule)
    - col_start: in {ctg_start,mol_start}
    - col_end:   in {ctg_end,mol_end}, assumed to be consistent with col_start
    Output:
    - list([start,end,ident]) where:
      - start, end are the 1-base coordinates of the interval
      - ident is the average identity over the interval defined as follows
        if the interval is composed of k mappings m1 to mk,
        ident = (sum_{i=1 to k} length(mi)*ident(mi) ) / (sum_{i=1 to k} length(mi))
        and the lngth of mi is taken either on the contig or te molecule as defined 
        by col_start,col_end.
    This implies that for assemblies for which contigs do overlap (such as SPAdes),
    bases will be considered several times.
    This is an approximation of the identity over the alignment defined by the
    mappings aggregated into an interval.
    '''
    ctg_mol_mappings_df.sort_values(by=col_end, inplace=True)
    ctg_mol_mappings_df.sort_values(by=col_start, inplace=True)
    # Resulting list of intervals, empty to start
    ctg_mol_intervals = []
    # Features of the current interval
    interval_start,interval_end = 0,-1
    interval_len,interval_ident = 1,0.0
    for _,row in ctg_mol_mappings_df.iterrows():
        mapping_len = row[col_end]-row[col_start]+1
        mapping_ident = mapping_len * row['ident']
        if row[col_start] > interval_end+1: # Gap
            # Closing and recording current interval
            ctg_mol_intervals.append([
                interval_start,interval_end,
                round(interval_ident/interval_len,4)
            ])
            # Starting a new current interval
            interval_start,interval_end = row[col_start],row[col_end]
            interval_len,interval_ident = mapping_len,mapping_ident
        else: # Overlap or contiguous mappings
            if row[col_end] > interval_end: 
                interval_end = row[col_end]
            interval_ident += mapping_ident
            interval_len += mapping_len
    # Recording current interval, not closed by starting a new one
    ctg_mol_intervals.append([
        interval_start,interval_end,
        round(interval_ident/interval_len,4)
    ])
    # Excluding first interval that was empty
    return ctg_mol_intervals[1:]
    
def compute_ctg_mol_pairs_intervals(in_mappings_df):
    '''
    Compute for each pair (mol,ctg) in in_mappings_df a list of triplets 
    (interval_start,interval_end,identity)
    corresponding of intervals along the contigs defined 
    by overlapping mappings, with identity averaged per base.
    Output:
    list (
      ctg_id,mol_id,mol_type,ctg_len,mol_len,
      list([ctg_start,ctg_end,ident]),list([mol_start,mol_end,ident])
    )
    '''
    grouped_df = in_mappings_df.groupby(
        by=['ctg_id','mol_id','mol_type','ctg_len','mol_len'], axis=0
    )
    ctg_mol_pairs_intervals = []
    for (ctg_id,mol_id,mol_type,ctg_len,mol_len),df_group in grouped_df:
        intervals_ctg_group = __compute_ctg_mol_intervals(
            df_group, 'ctg_start', 'ctg_end'
        )
        intervals_mol_group = __compute_ctg_mol_intervals(
            df_group, 'mol_start', 'mol_end'
        )
        ctg_mol_pairs_intervals.append([
            ctg_id,mol_id,mol_type,ctg_len,mol_len,
            intervals_ctg_group,intervals_mol_group
        ])
    return ctg_mol_pairs_intervals

INTERVALS_COL_NAMES = [
    'contig','molecule','molecule type','contig length',
    'molecule length','contig intervals','molecule intervals'
]
INTERVALS_COL_TYPES = {
    'contig': str, 'molecule': str, 'molecule type': str,
    'contig length': int,'molecule length': int,
    'contig intervals': str,'molecule intervals': str
}

def __intervals_to_str(intervals):
    return  INTERVALS_SEP.join(
        [INTERVAL_SEP.join([str(y) for y in x]) for x in intervals]
    )

def __intervals_from_str(intervals):
    return [
        [
            int(x.split(INTERVAL_SEP)[0]),
            int(x.split(INTERVAL_SEP)[1]),
            float(x.split(INTERVAL_SEP)[2])
        ]
        for x in intervals.split(INTERVALS_SEP)
    ]

def write_ctg_mol_pairs_intervals_TSV(intervals, out_tsv_file):
    '''
    Write contig/molecule pairs intervals into a TSV file
    '''
    formatted_intervals = [
        [
            x[0],x[1],x[2],x[3],x[4],
            __intervals_to_str(x[5]),
            __intervals_to_str(x[6])
        ]
        for x in intervals
    ]
    formatted_intervals.sort(key=lambda x: x[2])
    with open(out_tsv_file,'w') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(INTERVALS_COL_NAMES)
        writer.writerows(formatted_intervals)

def read_ctg_mol_pairs_intervals_TSV(in_tsv_file):
    '''
    Reads an intervals file and returns a dataframe with same columns
    where eac row corresponds to a single pair contig,molecule
    Intervals are split into a list of [start,end,identity]
    '''
    result = pd.read_csv(
            in_tsv_file, sep='\t', header=0, dtype=INTERVALS_COL_TYPES
        )        
    modify_column_function(
        result, 'contig intervals', lambda x: __intervals_from_str(x)
    )
    modify_column_function(
        result, 'molecule intervals', lambda x: __intervals_from_str(x)
    )
    return result

# Defining plASgraph ground truth

def __compute_ctg_cov_ident(
        in_ctg_mol_pair_row, min_interval_len, min_interval_ident
):
    '''
    in_ctg_mol_pair_row: intervals dataframe row corresponding to a single 
    pair contig,molecule
    min_interval_len, min_interval_ident: int,float
    computes coverage and identity from all contig intervals
    coverage: sum of lengths of intervals 
    (assumed to be non-overlapping)
    identity: average identity per base 
    '''
    def __test_interval(x):
        test_len = x[1]-x[0]+1 >= min_interval_len
        test_ident = x[2] >= min_interval_ident
        return test_len and test_ident
    ctg_intervals = [
        x
        for x in in_ctg_mol_pair_row['contig intervals']
        if __test_interval(x)
    ]
    if len(ctg_intervals) == 0: # All intervals discarded
        ctg_cov,ctg_ident=0,0.0
    else:
        ctg_cov = sum(
            [x[1]-x[0]+1 for x in ctg_intervals]
        )/in_ctg_mol_pair_row['contig length']
        ctg_ident = sum(
            [(x[1]-x[0]+1)*x[2] for x in ctg_intervals]
        )/ctg_cov
    return ctg_cov,ctg_ident

def compute_plasgraph_ground_truth(
        in_mappings_df, in_ctg_ids, 
        min_interval_len, min_interval_ident,
        min_ctg_cov, min_ctg_ident,
        out_intervals_file
):
    '''
    Computes a dict(ctg_id -> dict (mol_type: 0/1))
    for all contigs in ctg_ids using the intervals in 
    in_ctg_mol_pairs_intervals_df where each row corresponds
    to a single contig,molecule pair
    min_interval_len, min_interval_ident: int,float
    used to filter intervals that are too short or too low quality
    min_ctg_cov,min_ctg_ident: int,float
    thresholds to decide if a contig is positive (1) or negative (0)
    based on coverage and identity of all non-discarded intervals
    Intevals are written in TSV format in out_intervals_file
    '''
    ctg_score = {
        ctg_id: {MOL_CHR:0, MOL_PLS:0}
        for ctg_id in in_ctg_ids
    }
    ctg_mol_pairs_intervals = compute_ctg_mol_pairs_intervals(
        in_mappings_df
    )
    write_ctg_mol_pairs_intervals_TSV(
        ctg_mol_pairs_intervals, out_intervals_file
    )
    ctg_mol_pairs_intervals_df = read_ctg_mol_pairs_intervals_TSV(
        out_intervals_file
    )
    for _,ctg_mol_pair in ctg_mol_pairs_intervals_df.iterrows():
        ctg_id = ctg_mol_pair['contig']
        mol_type = ctg_mol_pair['molecule type']
        ctg_cov,ctg_ident = __compute_ctg_cov_ident(
            ctg_mol_pair, min_interval_len, min_interval_ident
        )
        if ctg_cov >= min_ctg_cov and ctg_ident >= min_ctg_ident:
            ctg_score[ctg_id][mol_type] = 1
    return ctg_score

PLASGRAPH_COL_NAMES = [
    'contig', 'label', 'length', 'chrom_score', 'plasmid_score'
]

def write_plasgraph_ground_truth(
        in_ctg_score_dict, in_ctg_len_dict, out_tsv_file
):
    '''
    in_ctg_scores_dict: dict((ctg_id -> dict (mol_type: 0/1))
    in_ctg_len_dict: dict(ctg_id -> ctg length)
    '''
    label = {
        (0,0): 'unlabeled', (1,0): 'chromosome',
        (0,1): 'plasmid', (1,1): 'ambiguous'
    }
    formatted_ctg_class = [
        [
            ctg_id, label[(ctg_score[MOL_CHR],ctg_score[MOL_PLS])],
            in_ctg_len_dict[ctg_id], ctg_score[MOL_CHR], ctg_score[MOL_PLS]
        ]
        for ctg_id,ctg_score in in_ctg_score_dict.items()
    ]
    with open(out_tsv_file,'w') as out_file:
        writer = csv.writer(out_file, delimiter=PLASGRAPH_SEP1)
        writer.writerow(PLASGRAPH_COL_NAMES)
        writer.writerows(formatted_ctg_class)

# Defining PlasBin ground truth

def compute_plasbin_ground_truth(
        in_mappings_df,
        min_interval_len, min_interval_ident,
        min_ctg_cov, min_ctg_ident
):
    '''
    Computes a list(
    [
      plasmid ID, contig ID, contig coverage, plasmid length, 
      contig length, intervals on molecule, intervals on contig
    ]
    for each pair plasmid,contig where from intervals 
    - of length >= min_interval_len and 
    - identity is >= min_interval_ident
    the contig coverage by intervals is >= min_ctg_cov and the has
    identity >= min_ctg_ident
    '''
    # Keeping only pairs contig,molecule where molecule is a plasmid
    
    ctg_mol_pairs_intervals = compute_ctg_mol_pairs_intervals(
        in_mappings_df
    )
    ctg_mol_pairs_intervals_df = pd.DataFrame(
        data=ctg_mol_pairs_intervals,
        columns=INTERVALS_COL_NAMES
    )
    plasbin_df = ctg_mol_pairs_intervals_df.loc[
        ctg_mol_pairs_intervals_df['molecule type']==MOL_PLS
    ]
    pls_ctg_list = []
    for _,row in plasbin_df.iterrows():
        ctg_cov,ctg_ident = __compute_ctg_cov_ident(
            row, min_interval_len, min_interval_ident
        )
        if ctg_cov >= min_ctg_cov and ctg_ident >= min_ctg_ident:
            pls_ctg_list.append([
                row['molecule'],row['contig'],ctg_cov,
                row['molecule length'],row['contig length'],
                __intervals_to_str(row['molecule intervals']),
                __intervals_to_str(row['contig intervals'])
            ])
    return pls_ctg_list

PLASBIN_COL_NAMES = [
    'plasmid', 'contig', 'contig coverage',
    'plasmid length', 'contig length',
    'plasmid intervals', 'contig intervals'
]

def write_plasbin_ground_truth(in_plasbin_list, out_tsv_file):
    with open(out_tsv_file,'w') as out_file:
        writer = csv.writer(out_file, delimiter=PLASBIN_SEP1)
        writer.writerow(PLASBIN_COL_NAMES)
        writer.writerows(in_plasbin_list)
