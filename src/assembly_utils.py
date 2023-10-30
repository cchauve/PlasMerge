#!/usr/bin/env python

'''
Manipulating FASTA and GFA files
'''

import csv
import gzip
import os
import argparse
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Misc

def __write_dict_TSV(in_dict, out_file_path):
    with open(out_file_path, 'w') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=list(in_dict[0].keys()), delimiter='\t')
        writer.writeheader()
        writer.writerows(in_dict)

## Reading FASTA files

def __read_FASTA(in_file_path, record_fun=lambda x: x):
    '''
    Returns a dictionary contig id -> record_fun(contig lengthrecord) (empty if file does not exist)
    TODO: use the io_input classes
    '''
    if os.path.isfile(in_file_path):
        return {
            id: record_fun(record)
            for id,record in SeqIO.to_dict(SeqIO.parse(gzip.open(in_file_path, 'rt'), 'fasta')).items()
        }
    else:
        return {}

def read_FASTA_id(in_file_path):
    '''
    Returns a list of records ID
    '''
    return list(__read_FASTA(in_file_path, record_fun=lambda x: None).keys())

def read_FASTA_len(in_file_path):
    '''
    Returns a dictionary contig id -> contig length (empty if file does not exist)
    '''
    return __read_FASTA(in_file_path, record_fun=lambda x: len(x.seq))

def read_FASTA_seq(in_file_path):
    '''
    Returns a dictionary contig id -> contig sequence (empty if file does not exist)
    '''
    return __read_FASTA(in_file_path, record_fun=lambda x: x.seq)

## Reading GFA files

def __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: x):
    '''
    Returns a dictionary vertex id -> ctg_fun(vertex)  (empty if file does not exist)
    TODO: use the io_input classes
    '''
    result = {}
    if os.path.isfile(in_file_path):
        with gzip.open(in_file_path, 'rt') as in_file:
            for line in in_file.readlines():
                line_split = line.strip().split('\t')
                if line[0] == 'S':
                    ctg_id,ctg_data = line_split[1],line_split[2:]
                    result[ctg_id] = ctg_fun(ctg_data)
    return result

def __GFA_S_attributes(in_data, in_key):
    '''
    Reads the list of attributes of a contig and returns the str value for attribute in_key, None if absent.
    '''
    attributes = {x.split(':')[0]: x.split(':')[2] for x in in_data}
    if in_key in attributes.keys():
        return attributes[in_key]
    else:
        return None

def read_GFA_ctg_id(in_file_path):
    '''
    Returns a list of contigs ID
    '''
    return list(__read_GFA_ctgs(in_file_path, ctg_fun=lambda x: None).keys())
    
def read_GFA_ctg_len(in_file_path):
    '''
    Returns a dictionary contig id -> contig length  (empty if file does not exist)
    '''
    return __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: len(x[0]))

def read_GFA_ctg_seq(in_file_path):
    '''
    Returns a dictionary contig id -> contig sequence  (empty if file does not exist)
    '''
    return __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: x[0])

def read_GFA_ctg_RC(in_file_path):
    '''
    Returns a dictionary contig id -> contig read count  (empty if file does not exist)
    '''
    return __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: int(__GFA_S_attributes(x, 'RC')))

def read_GFA_ctg_KC(in_file_path):
    '''
    Returns a dictionary contig id -> contig kmer count  (empty if file does not exist)
    '''
    return __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: int(__GFA_S_attributes(x, 'KC')))

def read_GFA_ctg_LN(in_file_path):
    '''
    Returns a dictionary contig id -> contig read countlen  (empty if file does not exist)
    '''
    return __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: int(__GFA_S_attributes(x, 'LN')))

def read_GFA_ctg_dp(in_file_path):
    '''
    Returns a dictionary contig id -> contig read countlen  (empty if file does not exist)
    '''
    return __read_GFA_ctgs(in_file_path, ctg_fun=lambda x: float(__GFA_S_attributes(x, 'dp')))

def __read_GFA_edges(in_file_path):
    '''
    Returns a dictionary ctg1_id -> list (ctg2_id, sign1, sign2) describing all edges
    (empty if no file)
    '''
    result = {}
    if os.path.isfile(in_file_path):
        with gzip.open(in_file_path, 'rt') as in_file:
            for line in in_file.readlines():
                line_split = line.strip().split('\t')
                if line[0] == 'L':
                    ctg1_id,sign1,ctg2_id,sign2 = line_split[1:5]
                    result[ctg1_id] = (ctg2_id,sign1,sign2)
    return result

def write_GFA_to_FASTA(in_GFA_file, out_FASTA_file):
    '''
    Write the contigs of a GFA file in FASTA format
    '''
    GFA_ctg_seqs = read_GFA_ctg_seq(in_GFA_file)
    ctg_records = [
        SeqRecord(Seq(y), id=x, name=x, description=f'{x}.GFA')
        for x ,y in GFA_ctg_seqs.items()
    ]
    with gzip.open(out_FASTA_file, 'wt') as out_file:
        SeqIO.write(ctg_records, out_file, 'fasta')

def extract_matches_from_BLAST(
        in_FASTA_file, in_GFA_file, in_BLAST_file,
        in_format=6, min_cov=0.95, min_ident=0.95
):
    '''
    From the result of blast over in_FASTA_file, in_GFA_file 
    identify for every contig in in_FASTA_file 
    all contigs in in_GFA_file such that the alignment covers
    min_cov of both sequences, with min_ident
    in_format: BLAST file format
    return a dict FASTA_ctg_id -> list(GFA_ctg_id,FASTA_len,GFA_len,alignment_len,alignmend_identity)
    '''
    FASTA_ctg_len = read_FASTA_len(in_FASTA_file)
    GFA_ctg_len = read_GFA_ctg_len(in_GFA_file)
    with open(in_BLAST_file) as in_file:
        result = {ctg_id: [] for ctg_id in FASTA_ctg_len.keys()}
        for blast_result in in_file.readlines():
            blast_record = blast_result.rstrip().split('\t')
            qseqid,sseqid = blast_record[0],blast_record[1]
            pident = float(blast_record[2])/100.0
            sqlen = float(blast_record[3])
            qcov = sqlen/FASTA_ctg_len[qseqid]
            scov = sqlen/GFA_ctg_len[sseqid]
            if pident >= min_ident and qcov >= min_cov and scov >= min_cov:
                result[qseqid].append((sseqid,FASTA_ctg_len[qseqid],GFA_ctg_len[sseqid],sqlen,pident))
    return result
    
def export_matches_from_BLAST(in_matches_dict, out_TSV_file):
    '''
    Creates a TSV file FASTA ctg id<TAB>comma-separated list of matched GFA ctg id/NA if none
    '''
    out_dict = [
        {
            'FASTA_ctg':
            x,
            'GFA_ctg:FASTA_len:GFA_len:alg_len:indentity':
            'NA' if len(z)==0 else ','.join([':'.join([str(u) for u in y]) for y in z])
        }
        for x,z in in_matches_dict.items()
    ]    
    __write_dict_TSV(out_dict, out_TSV_file)    

def main():
    parser = argparse.ArgumentParser(description='Manipulating assemblies.')
    subparsers = parser.add_subparsers(help='sub-command help')
    # convert command arguments
    convert_parser = subparsers.add_parser('convert')
    convert_parser.set_defaults(cmd='convert')
    convert_parser.add_argument('input', type=str, help='Input gzipped GFA file')
    convert_parser.add_argument('output', type=str, help='Output gzipped FASTA file')
    # match command arguments
    match_parser = subparsers.add_parser('match')
    match_parser.set_defaults(cmd='match')
    match_parser.add_argument('input_FASTA', type=str, help='Input gzipped FASTA file')
    match_parser.add_argument('input_GFA', type=str, help='Input gzipped GFA file')
    match_parser.add_argument('input_BLAST', type=str, help='Input BLAST results file')
    match_parser.add_argument('output', type=str, help='Output matches TSV file')
    match_parser.add_argument('min_cov', type=float, help='Minimum coverage')
    match_parser.add_argument('min_ident', type=float, help='Minimum identity')    

    args = parser.parse_args()

    if args.cmd == 'convert':
        write_GFA_to_FASTA(args.input, args.output)
    elif args.cmd == 'match':
        matches_dict = extract_matches_from_BLAST(
            args.input_FASTA, args.input_GFA, args.input_BLAST,
            in_format=6, min_cov=args.min_cov, min_ident=args.min_ident
        )
        export_matches_from_BLAST(matches_dict, args.output)

if __name__ == "__main__":
    main()
