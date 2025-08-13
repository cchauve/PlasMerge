#!/usr/bin/python

# TO DO
# Rewrite PBF_output.parse_bins_file

import pandas as pd
import networkx as nx

from gfa_fasta_utils import (
    read_GFA_id,
    read_GFA_seq,
    read_GFA_len,
    read_GFA_normalized_coverage,
    read_GFA_links,
    GFA_FROM_ORIENT_KEY,
    GFA_TO_KEY,
    GFA_TO_ORIENT_KEY
)

from PBF_utils import (
    read_gc_intervals_file,
    compute_gc_probabilities,
    read_pls_score_file,
    GC_COUNT_KEY,
    GC_RATIO_KEY,
    LENGTH_KEY,
    DEFAULT_GC_INTERVALS
)

# Assembler used to generate the GFA file
UNICYCLER_TAG = 'unicycler'
SKESA_TAG = 'skesa'

"""
Class recording a contig and its attributes:
- sequence
- sequence length
- read depth
- GC content
"""
class Contig():
    def __init__(self, seq, length, rd):
        """
        Initiate a Contig object
        Args:
           - seq (str): contig sequence
           - len (int): contig length
           - rd (float): rad depth
        Returns:
           NA
        """
        
        self.seq = seq
        self.len = length
        self.rd = rd
        self.gc = (seq.count('G') + seq.count('C'))
        
    def set_rd(self, rd):
        """
        Stores the read depth of the contig sequence (float)
        """
        self.rd = float(rd)	

    def get_seq(self):
        """
        Returns the sequence of the contig (str)
        """
        return self.seq

    def get_len(self):
        """
        Returns the length of the contig (int)
        """		
        return self.len

    def get_rd(self):
        """
        Returns the read depth of the contig (float)
        """		
        return self.rd

    def get_gc(self):
        """
        Returns the GC content of the contig (int)
        """		
        return self.gc

    def get_gc_ratio(self):
        """
        Returns the GC content ratio of the contig (float in [0,1])
        """		
        return (self.gc / self.len if self.len > 0 else 0.0)
    
# TAG of GFA file for recording normalized coverage
ASSEMBLER_COV_TAG = {
    UNICYCLER_TAG: 'dp',
    SKESA_TAG: None
}
# Contigs extremities conversion
CTG_FROM_EXT = {'+': 'h', '-': 't'}
CTG_TO_EXT = {'+': 't', '-': 'h'}    

"""
Class recording an assembly graph encoded as 
- a dictionary of Contig objects
  Dictionary(contig id -> Contig object)
- a list of directed edges
  List((contig id (str),{'h','t'}),(contig id (str),{'h','t'}))
  Contig extremity from, contig extremity to (edges are directed)
"""
class Assembly():
    def __init__(self, assembly_graph_file, gzipped=False, assembler=UNICYCLER_TAG):
        """
        Initiate an Assembly object from a GFA file
        Args:
        - assembly_graph_file (str): path to a GFA file
        - gzipped (bool): True if GFA file is gzippd
        - assembler (str): UNICYCLER_TAG or SKESA_TAG
        Returns:
        NA
        """
        self.ctg_dict = {}
        self.edge_list = []
        self.parse_gfa_file(assembly_graph_file, gzipped, assembler)

    def parse_gfa_file(self, assembly_graph_file, gzipped, assembler):
        """
        Populate attributes from a GFA file
        Args:
        - assembly_graph_file (str): path to a GFA file
        - gzipped (bool): True if GFA file is gzippd
        - assembler (str): UNICYCLER_TAG or SKESA_TAG
        Returns:
        NA
        """
        # Reading contigs
        ctgs_id = read_GFA_id(assembly_graph_file, gzipped=gzipped)
        ctgs_seq = read_GFA_seq(assembly_graph_file, gzipped=gzipped)
        ctgs_len = read_GFA_len(assembly_graph_file, gzipped=gzipped)
        ctgs_rd = read_GFA_normalized_coverage(
            assembly_graph_file,
            cov_key=ASSEMBLER_COV_TAG[assembler],
            gzipped=gzipped
        )
        for ctg_id in ctgs_id:
            self.ctg_dict[ctg_id] = Contig(
                ctgs_seq[ctg_id], ctgs_len[ctg_id], ctgs_rd[ctg_id]
            )
        # Reading edges
        edges_dict = read_GFA_links(assembly_graph_file, gzipped=gzipped)
        for ctg_id,ctg_edges_from_list in edges_dict.items():
            for edge in ctg_edges_from_list:
                edge_from = (ctg_id, CTG_FROM_EXT[edge[GFA_FROM_ORIENT_KEY]])
                edge_to = (edge[GFA_TO_KEY], CTG_TO_EXT[edge[GFA_TO_ORIENT_KEY]])
                self.edge_list.append((edge_from, edge_to))

    def get_ctg(self, ctg_id):
        """
        Returns Contig object for contig ctg_id (str)
        """
        return self.ctg_dict[ctg_id]

    def get_ctg_ids(self):
        """
        Returns list of contigs ID (str)
        """
        return list(self.ctg_dict.keys())
	 
    def get_edges(self):
        """
        Returns list of edges:
        List((contig id (str),{'h','t'}),(contig id (str),{'h','t'}))
        """
        return self.edge_list
	
    def get_nctgs(self):
        """
        Returns number of contigs in the graph (int)
        """
        return len(self.ctg_dict.keys())
	
    def get_nedges(self):
        """
        Returns number of edges in the graph (int)
        """
        return len(self.edge_list)

    def get_ctg_lens(self):
        """
        Returns contig lengths:
        Dictionary(contig id (str) -> contig length (int))
        """
        return {
            ctg_id: ctg.get_len()
            for ctg_id,ctg in self.ctg_dict.items()
        }

    def get_ctg_rds(self):
        """
        Returns contig lengths:
        Dictionary(contig id (str) -> contig read depth (float))
        """
        return {
            ctg_id: ctg.get_rd()
            for ctg_id,ctg in self.ctg_dict.items()
        }    

    def to_graph(self, directed=False):
        """
        Returns a networkx graph encoding the assembly graph
        """
        G = nx.Graph()
        ctg_list = self.get_ctg_ids()
        G.add_nodes_from(ctg_list)
        edge_list = [(edge[0][0], edge[1][0]) for edge in self.get_edges()]
        G.add_edges_from(edge_list)
        if directed:
            return G.to_directed()
        else:
            return G

def _read_gc_intervals(gc_intervals_file):
    """
    Reads a GC intervals file.
    Args:
    gc_intervals_file:
    - (str) : path to a GC intervals file
    - or None
    GC intervals file format: increasing list of numbers from 0 to 1.
    Warning: 
    file format is not checked
    Returns:
    List(float): increasing list of numbers from 0 to 1
    Represents boundaries of GC ratio intervals.
    """
    if gc_intervals_file is not None:
        return read_gc_intervals_file(gc_intervals_file)
    else:
        return DEFAULT_GC_INTERVALS

"""
Class recording the features of a set of plasmid bins
- pls_bins: contig content of each bin
  Dictionary(plasmid id (str) -> plasmid content, List((contig id (str), contig multiplicity (int)))
  default multiplicity: 1
- copy_number: copy number of each plasmid (0.0 if not specified)
  Dictionary(plasmid id (str) -> float)
- source: source of the plasmid bins (ground truth or prediction tool)
  (str)
"""
class Plasmids:
    def __init__(self, source, pls_file, assembly=None):
        """
        Initiate a Plasmids object
        Args:
        - pls_file: (str) plasmid bins file
        - assembly: (Assembly) assembly object used to associate a copy number if not given
        Format: see parse_plasmids_file below
        source:     (str), source of plasmid / plasmid bins, either ground truth or name of tool
        """
        self.parse_plasmids_file(pls_file, assembly=assembly)
        self.source = source

    def parse_plasmids_file(self, pls_file, assembly=None):
        """
        Reading a plasmids bins file.
        Format:
        - line 1: header "plasmid contigs <optional:copy_number>"
        - plasmid line: <plasmid id>TAB<comma-spearated list of contig_id:multiplicity>TAB<copy number>
        Args: 
        - pls_file: (str) path to plasmid file
        - assembly: (Assembly) assembly object used to associate a copy number if not given
        Returns:
        NA
        """
        plasmids_df = pd.read_csv(
            pls_file, sep='\t', header=0, skip_blank_lines=True
        )
        self.pls_bins = {}
        self.copy_number = {}
        for idx,row in plasmids_df.iterrows():
            pls_id = row['plasmid']
            ctgs = row['contigs'].split(',')
            self.pls_bins[pls_id] = [
                (ctg.rsplit(':',1)[0],int(float(ctg.rsplit(':',1)[1])))
                for ctg in ctgs
            ]
            if assembly is not None:
                # Min read depth of all contigs in the plasmid bin
                self.copy_number[pls_id] = min([
                    assembly.get_ctg(ctg_id).get_rd()
                    for (ctg_id,_) in self.pls_bins[pls_id]
                ])
            elif 'copy_number' in row.keys():
                self.copy_number[pls_id] = float(row['copy_number'])
            else:
                self.copy_number[pls_id] = 0.0

    def get_pls_ids(self):
        """
        Return list of predicted plasmids ID (List(str))
        """
        return list(self.pls_bins.keys())
    
    def get_pls_content(self, with_mult=True):
        """
        Return bin contig content for all plasmids
        Dictionary(plasmid ID (str) ->
        if with_mult is True: (List(contig id (str), multiplicity (int)))
        if with_mult is False: (List(contig id (str))))
        """
        if with_mult:
            return self.pls_bins
        else:
            return {
                pls_id: [ctg[0] for ctg in self.pls_bins[pls_id]]
                for pls_id in self.pls_bins.keys()
            }

    def get_pls_copy_number(self):
        """
        Return
        Dictionary(plasmid id (str) -> float)
        """
        return self.copy_number
        
    def get_source(self):
        """
        Return plasmids source (str)
        """
        return self.source

def _compute_optimal_gc_bin(pls_bins, gc_probs, nb_gc_intervals):
    """
    Compute the optimal GC interval for a set of plasmid bins
    Args
    - pls_bins: (Plasmids) plasmid bins
    - gc_probs: Dictionary(contig id (str) -> List(float)) GC proba of each contig/GC bin of input
    - nb_gc_intervals: (int) number of GC intervals
    Returns
    - Dictionary(plasmid id (str) -> (int) index of optimal GC interval
    """
    # Computing GC penalty for each contig and GC bin as defined in PlasBin-flow
    ctg_gc_penalties = {}
    for ctg_id,ctg_gc_probs in gc_probs.items():
        ctg_gc_prob_max = max(ctg_gc_probs)
        ctg_gc_penalties[ctg_id] = [
            - (ctg_gc_prob_max - ctg_gc_prob)
            for ctg_gc_prob in ctg_gc_probs
        ]
    # Computing the GC bin with optimal penalty
    pls_opt_gc_bin = {}
    plasmids = pls_bins.get_pls_content(with_mult=True)
    for pls_id,pls_bin in plasmids.items():
        pls_opt_penalty = -len(pls_bin)*1.0
        for i in range(nb_gc_intervals):
            pls_gc_penalty = sum([
                ctg_mult*ctg_gc_penalties[ctg_id][i]
                for (ctg_id,ctg_mult) in pls_bin
            ])
            if pls_gc_penalty > pls_opt_penalty:
                pls_opt_penalty = pls_gc_penalty
                pls_opt_gc_bin[pls_id] = i
    return pls_opt_gc_bin

"""
Class recording the input to PlasMerge
assembly:     Assembly object (assembly graph)
assembler:    UNICYCLER_TAG or SKESA_TAG         
gc_intervals: List(float) list of GC interval boundaries
gc_probs:     Dictionary(contig id (str) -> GC proba for each interval List(float))
pls_score:    Dictionary(contig id (str) -> plasmid score (float))
pls_bins:     Plasmids object (plasmid bins)
gc_bins:      (int) index of the optimal GC bin per plasmid bin
""" 
class PBM_input():
    def __init__(self,
                 assembly_graph_file,
                 pls_score_file,
                 gc_intervals_file,
                 pls_bins_file,
                 source,
                 gzipped=False,
                 assembler=UNICYCLER_TAG
    ):
        """
        Initiate PBF_input object
        Args:
        - assembly_graph_file (str): path to a GFA file
        - pls_score_file (str):      path to a plasmid score file (format: ctg_id<TAB>score)
        - gc_intervals_file (str):   path to GC intervals file (format: see _read_gc_intervals)
        - pls_bins_file (str):       path to a plasmid bins file (format: see Plasmids)
        - source: (str)              source of the plasmid bins (ground truth or prediction tool)
        - gzipped (bool):            True if GFA file is gzippd
        - assembler (str):           UNICYCLER_TAG or SKESA_TAG
        Returns:
        NA
        """
        self.assembly = Assembly(assembly_graph_file, gzipped=gzipped, assembler=assembler)
        self.pls_score = read_pls_score_file(pls_score_file)
        self.gc_intervals = _read_gc_intervals(gc_intervals_file)
        nb_gc_intervals = len(self.gc_intervals) - 1
        ctgs_gc = {}
        for ctg_id in self.assembly.get_ctg_ids():
            ctg = self.assembly.get_ctg(ctg_id)
            ctgs_gc[ctg_id] = {
                GC_COUNT_KEY: ctg.get_gc(),
                GC_RATIO_KEY: ctg.get_gc_ratio(),
                LENGTH_KEY: ctg.get_len()
            }
        self.gc_probs = compute_gc_probabilities(ctgs_gc, self.gc_intervals)
        self.pls_bins = Plasmids(source, pls_bins_file, assembly=self.assembly)
        self.gc_bins = _compute_optimal_gc_bin(self.pls_bins, self.gc_probs, nb_gc_intervals)
        
    def get_gc_intervals(self):
        """
        Returns list of GC interval boundaries:
        List(float) sorted increasingly
        """
        return self.gc_intervals

    def get_gc_probs(self):
        """
        Returns GC probabilities:
        Dictionary(contig id (str) -> List(float))
        where element in position i in the list is the probability to be in the (i+1)th interval
        """
        return self.gc_probs

    def get_pls_scores(self):
        """
        Returns plasmid scores:
        Dictionary(contig id (str)-> plasmid score (float))
        """
        return self.pls_score

    def get_assembly(self):
        """
        Returns assembly graph (Assembly object)
        """
        return self.assembly

    def get_pls_bins(self):
        """
        Returns plasmid bins (Plasmids object)
        """
        return self.pls_bins

    def get_gc_bins(self):
        """
        Returns Dictionary(plasmid id (str) -> int (index of optimal GC bin)
        """
        return self.gc_bins

# # Testing the classes above
# if __name__ == "__main__":
#     import os
    
#     samples = ['SAMN32247302', 'SAMN32247345', 'SAMN32247425', 'SAMN32247519', 'SAMN32247522']
#     root = os.path.normpath('../test')
#     gc_intervals_file = os.path.join(root, 'gc_intervals.txt')
#     sources = ['gt', 'mob', 'gp', 'pbf']
    
#     for sample in samples:
#         print(f'SAMPLE: {sample}')
    
#         gfa_file = os.path.join(root, 'gfas', f'{sample}.assembly.gfa.gz')
#         pls_scores_file = os.path.join(root, 'scores', f'{sample}.scores.tsv')
        
#         assembly = Assembly(gfa_file, gzipped=True)
#         print(f'\tGFA\tnb contigs\t{assembly.get_nctgs()}')
#         print(f'\tGFA\tnb edges\t{assembly.get_nedges()}')
#         ctg_id = assembly.get_ctg_ids()[0]
#         ctg = assembly.get_ctg(ctg_id)
#         print(f'\tCONTIG {ctg_id} {ctg.get_len()} {ctg.get_rd()} {ctg.get_gc()} {ctg.get_gc_ratio()}')
#         print(f'\tConversion to networkx graph')
#         G = assembly.to_graph()

#         gc_intervals = _read_gc_intervals(gc_intervals_file)
#         nb_gc_intervals = len(gc_intervals) - 1

#         for source in sources:
#             print(f'\tSOURCE: {source}')
            
#             pls_bins_file = os.path.join(root, 'pls_bins', f'{sample}.{source}.tsv')

#             pbm_input = PBM_input(gfa_file, pls_scores_file, gc_intervals_file, pls_bins_file, source, gzipped=True)            

#             if source == 'gt':
#                 print(f'\tGC INTERVALS\t{pbm_input.get_gc_intervals()}')
#                 gc_probs = pbm_input.get_gc_probs()
#                 scores = pbm_input.get_pls_scores()
#                 assembly = pbm_input.get_assembly()
#                 ctg_id = assembly.get_ctg_ids()[0]
#                 print(f'\tGC PROBA CONTIG {ctg_id}\t{gc_probs[ctg_id]}')
#                 print(f'\tSCORE CONTIG {ctg_id}\t {scores[ctg_id]}')
#             pls_bins = pbm_input.get_pls_bins()
#             print(f'\tSOURCE\t{pls_bins.get_source()}')
#             print(f'\tPLASMID CONTENT\t{pls_bins.get_pls_content()}')
#             print(f'\tPLASMID COPY NUMBER\t{pls_bins.get_pls_copy_number()}')            
#             print(f'\tOPT GC BINS\t{pbm_input.get_gc_bins()}')

