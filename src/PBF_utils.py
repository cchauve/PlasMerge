""" 
Functions and variables imported from PlasBin-flow 1.0.2, with logging removed
"""

import os
import scipy.special as sc
import scipy.integrate as integrate
import math

# Keys of dictionary associated to a FASTA record
GC_COUNT_KEY = 'gc_count'
GC_RATIO_KEY = 'gc_ratio'
LENGTH_KEY = 'length'

# Default GC content ratio intervals
DEFAULT_GC_INTERVALS = [0, 0.4,  0.45, 0.5, 0.55, 0.6, 1]

# Default plasmid score offset to compute the objective function
DEFAULT_SCORE_OFFSET = 0.5

def read_gc_intervals_file(gc_intervals_file):
    """
    Read GC intervals file
    Args:
        - gc_intervals_file (str): path to GC intervals file
    Returns:
        List(float): sorted list of GC intervals boundaries
    """
    with open(gc_intervals_file) as in_file:
        intervals = [
            float(x.rstrip())
            for x in in_file.readlines()
        ]
    return intervals

def compute_gc_probabilities(contigs_gc, gc_intervals):
    """
    Computes GC probabilities for a set of contigs and GC intervals
    Args:
        - contigs_gc (Dictionary): contig id -> {
          GC_COUNT_KEY: GC count, GC_RATIO_KEY: GC content ratio, LENGTH_KEY: contig length
          }
        - gc_intervals (List(float)): GC intervals boundaries

    Returns:
        (Dictionary) contig id -> List(float)
        where element in position i in the list is the probability to be in the (i+1)th interval
    """
    def _gprob2(n, g, p, m):
        """
        Compute probability of observing g GC nucleotides in a contig
        of length n within a molecule of GC content p using pseducount m.
        Done via logarithm to avoid overflow.
        """
        def _combln(n, k):
            """ Compute ln of n choose k. Note than n! = gamma(n+1) """
            return sc.gammaln(n + 1) - sc.gammaln(k + 1) - sc.gammaln(n - k + 1)
        alpha = m*p
        beta = m*(1-p)
        resultln = _combln(n, g) + sc.betaln(g + alpha, n - g + beta) - sc.betaln(alpha, beta)
        return math.exp(resultln)
    
    m = 10
    ctgs_gcp = {}
    for ctg_id,ctg_data in contigs_gc.items():
        n = ctg_data[LENGTH_KEY]
        g = ctg_data[GC_COUNT_KEY]
        total = 0
        gcp_array = []
        for i in range(0, len(gc_intervals)-1):
            gp2 = integrate.quad(
                lambda x: _gprob2(n,g,x,m),
                    gc_intervals[i], gc_intervals[i+1]
            )
            gp2 = gp2[0]/(gc_intervals[i+1] - gc_intervals[i])
            total += gp2
            gcp_array.append(gp2)
        ctgs_gcp[ctg_id] = [gcp/total for gcp in gcp_array]
    return ctgs_gcp

def read_pls_score_file(in_pls_score_file):
    """
    Reads a plasmid score file
    Args:
        in_pls_score_file (str): path to a plasmid score file
    Returns:
        Dictionary contig (str): plasmid score (float)
    Assumption:
        File existence and non-emptyness has been checked
    """
    pls_scores_dict = {}
    with open(in_pls_score_file) as in_file:
        for score_line in in_file.readlines():
            line = score_line.rstrip()
            line_split = line.split('\t')
            ctg_id = line_split[0]
            score = float(line_split[1])
            pls_scores_dict[ctg_id] = score
    return pls_scores_dict

