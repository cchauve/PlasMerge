"""
Computing the optimal merging for all pairs of bins for a given sample
"""

import os
from itertools import combinations
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import pandas as pd

from io_utils import PBM_input
from PBF_utils import DEFAULT_SCORE_OFFSET

# sums up GC and gene density terms for contigs in bin_ctgs
# into single node with name node
def collapse_obj_terms(pls_score, gc_probs, gc_ints, bin_ctgs, node):
    
    merged_gc = {}
    merged_pls_score = {}
    
    for ctg in bin_ctgs:
        max_prob = max(gc_probs[ctg])
        if merged_pls_score.get(node) != None:
            for i in gc_ints:
                merged_gc[(node, i)] += gc_probs[ctg][i] - max_prob
            merged_pls_score[node] += pls_score[ctg] - DEFAULT_SCORE_OFFSET
        else:
            for i in gc_ints:
                merged_gc[(node, i)] = gc_probs[ctg][i] - max_prob
            merged_pls_score[node] = pls_score[ctg] - DEFAULT_SCORE_OFFSET

    return merged_gc, merged_pls_score
        
# computes and returns objective function terms for each node 
# in modified graph (after contracting both bins into one node)
# to be used in post-processing MILP
def obj_terms(pls_score, gc_probs, gc_ints, bin1_ctgs, bin2_ctgs):
    
    bin1_gc, bin1_pls_score = collapse_obj_terms(pls_score, gc_probs, gc_ints, bin1_ctgs, 's')
    bin2_gc, bin2_pls_score = collapse_obj_terms(pls_score, gc_probs, gc_ints, bin2_ctgs, 't')
    bin1_gc.update(bin2_gc)
    bin1_pls_score.update(bin2_pls_score)
    
    for ctg in set(pls_score.keys()).difference(set(bin1_ctgs) | set(bin2_ctgs)):
        max_prob = max(gc_probs[ctg])
        for i in gc_ints:
            bin1_gc[(ctg, i)] = gc_probs[ctg][i] - max_prob
        bin1_pls_score[ctg] = pls_score[ctg] - DEFAULT_SCORE_OFFSET
      
    return bin1_gc, bin1_pls_score

# Removes read_depth used up by other bins
def remove_plasmid_rd(rd, flows, pls_bins, pls_ids, bin1, bin2):

    for key in pls_ids:
        if key in [bin1, bin2]:
            continue
        bin_mults = dict(pls_bins[key])
        for ctg in bin_mults:
            if rd.get(ctg, None) == None:
                continue
            else:
                rd[ctg] = rd[ctg] - (bin_mults[ctg] * flows[key])
                
# Removing vertices that now have read-depth less than threshold                
def remove_low_rd(gc_pen, pls_score_pen, rd, G, gc_ints, threshold):
    
    to_remove = [node for node in rd if rd[node] < threshold]
    for node in to_remove:
        for i in gc_ints:
            gc_pen.pop((node, i))
        pls_score_pen.pop(node)
        rd.pop(node)
        G.remove_node(node)

# contracting bins down to single vertices in networkx graph
def contracted_bins(graph, bin1_ctgs, bin2_ctgs):

    H = nx.DiGraph(graph)
    bin1_set = set(bin1_ctgs)
    bin2_set = set(bin2_ctgs)
    shared_ctgs = list(bin1_set & bin2_set)
    bin1_diff = list(bin1_set - bin2_set)
    bin2_diff = list(bin2_set - bin1_set)

    for ctg in bin1_diff[1:]:
        nx.contracted_nodes(H, bin1_diff[0], ctg, self_loops=False, copy=False)
    for ctg in bin2_diff[1:]:
        nx.contracted_nodes(H, bin2_diff[0], ctg, self_loops=False, copy=False)
    H = nx.relabel_nodes(H, {bin1_diff[0]:'s', bin2_diff[0]:'t'})
    if shared_ctgs:
        for ctg in shared_ctgs[1:]:
            nx.contracted_nodes(H, shared_ctgs[0], ctg, self_loops=False, copy=False)
        shared_edges = nx.edge_boundary(H, [shared_ctgs[0]], set(H.nodes()) - {'s', 't'})
        for source, nbr in shared_edges:
            H.add_edges_from([('s', nbr), (nbr, 's'), ('t', nbr), (nbr, 't')])
        H.add_edges_from([('s', 't'), ('t', 's')])
        H.remove_node(shared_ctgs[0])

    return H

# Creates MILP gurobipy object to be optimized when merging two bins
def generate_model(H, gc_pen, pls_score_pen, rd, gc_ints):
    
    model = gp.Model("model")
    model.setAttr('ModelSense', GRB.MAXIMIZE)
    edge_list = list(H.edges())
    node_list = list(H.nodes())

    x = model.addVars(edge_list, vtype=GRB.BINARY, name="x")
    q = model.addVars(edge_list, lb=0, vtype=GRB.CONTINUOUS, name="q")
    z = model.addVars(node_list, obj=pls_score_pen, vtype=GRB.BINARY, name="z")
    d = model.addVar(obj=1.0, lb=0, ub=10000, vtype=GRB.CONTINUOUS, name='d')
    GC = model.addVars(gc_ints, vtype=GRB.BINARY, name='GC')
    ctg_GC = model.addVars(z, GC, obj=gc_pen, vtype=GRB.BINARY, name='ctg_GC')

    # Dummy variable constraints so that d acts as the minimum read depth
    model.addConstrs(d <= 10000 - (10000 - rd[ctg]) * z[ctg] for ctg in node_list)
    # constraints to handle product of binary variables as one
    model.addConstrs(ctg_GC[ctg, i] <= z[ctg] for ctg in node_list for i in gc_ints)
    model.addConstrs(ctg_GC[ctg, i] <= GC[i] for ctg in node_list for i in gc_ints)
    model.addConstrs(ctg_GC[ctg, i] >= z[ctg] + GC[i] - 1 for ctg in node_list for i in gc_ints)
    model.addConstr(gp.quicksum(GC[i] for i in gc_ints) == 1)
    # Generic ESPP constraints
    model.addConstr(gp.quicksum(x['s', j] for j in H.successors('s')) - gp.quicksum(x[i, 's'] for i in H.predecessors('s')) == 1)
    model.addConstr(gp.quicksum(x['t', j] for j in H.successors('t')) - gp.quicksum(x[i, 't'] for i in H.predecessors('t')) == -1)
    model.addConstrs(gp.quicksum(x[k, j] for j in H.successors(k)) - gp.quicksum(x[i, k] for i in H.predecessors(k)) == 0 \
                     for k in node_list if (k != 's' and k != 't'))
    model.addConstrs(gp.quicksum(x[i, j] for j in H.successors(i)) <= 1 for i in node_list)
    # auxiliary subtour constraints
    model.addConstrs(q[i, j] <= (nx.number_of_nodes(H) - 1) * x[i, j] for (i, j) in edge_list)
    model.addConstr((gp.quicksum(q['s', j] for j in H.successors('s'))) == (gp.quicksum(z[i] for i in node_list if i != 's')))
    model.addConstrs(gp.quicksum(q[i, k] for i in H.predecessors(k)) - gp.quicksum(q[k, j] for j in H.successors(k)) == z[k] \
                     for k in node_list if k != 's')
    model.addConstrs(gp.quicksum(x[i, k] for i in H.predecessors(k)) == z[k] for k in node_list if k != 's')
    model.addConstr(z['s'] == 1)
    model.addConstr(z['t'] == 1)
    model.addConstr(gp.quicksum(q[i, 't'] for i in H.predecessors('t')) == 1)
    
    return model
    
# Outputs gurobipy model object corresponding to a post-processing
# MILP instance for the given pair of bins
def merging_model(pbm_input, bin1, bin2, threshold):
    # Reading the input files
    assembly = pbm_input.get_assembly()
    G = assembly.to_graph()
    pls_score = pbm_input.get_pls_scores()
    gc_probs = pbm_input.get_gc_probs()
    pls_bins = pbm_input.get_pls_bins()
    gc_ints = range(len(pbm_input.get_gc_intervals()) - 1)

    # Reading the set of all plasmid bins
    pls_ids = pls_bins.get_pls_ids()
    pls_content = pls_bins.get_pls_content(with_mult=True)
    flows = pls_bins.get_pls_copy_number()

    # Reading the two considered bins
    bin1_ctgs = [ctg_id for (ctg_id,_) in pls_content[bin1]]
    bin2_ctgs = [ctg_id for (ctg_id,_) in pls_content[bin2]]
    rd = {}

    # Processing the two bins to create a model
    gc_penalty, pls_score_penalty = obj_terms(pls_score, gc_probs, gc_ints, bin1_ctgs, bin2_ctgs)

    # Updating read depth
    for ctg_id in set(pls_score.keys()).difference(set(bin1_ctgs) | set(bin2_ctgs)):
        rd[ctg_id] = assembly.get_ctg(ctg_id).get_rd()
    rd['s'] = float(flows[bin1])
    rd['t'] = float(flows[bin2])
    remove_plasmid_rd(rd, flows, pls_content, pls_ids, bin1, bin2)
    H = contracted_bins(G, bin1_ctgs, bin2_ctgs)
    if threshold == 'min':
        remove_low_rd(gc_penalty, pls_score_penalty, rd, H, gc_ints, min(rd['s'], rd['t']))
    else:
        remove_low_rd(gc_penalty, pls_score_penalty, rd, H, gc_ints, threshold)
    
    return generate_model(H, gc_penalty, pls_score_penalty, rd, gc_ints)

GUROBI_ZERO = 0.0001

OUT_COLUMNS = [
    'SAMPLE', 'BIN1', 'BIN1_CONTIGS', 'BIN1_GC_BIN', 'BIN1_GC', 'BIN1_GD', 'BIN1_RD',
    'BIN2', 'BIN2_CONTIGS', 'BIN2_GC_BIN', 'BIN2_GC', 'BIN2_GD', 'BIN2_RD', 'MERGE_GC_BIN',
    'MERGE_GC', 'MERGE_GD', 'MERGE_RD', 'MERGE_PATH', 'MILP_INFEASIBLE', 'BINS_OVERLAP'
]

# returns relevant data from post-processing results
def merging_results(model, pbm_input, bin1, bin2):

    # returns string containing contig ID
    # and given properties of that contig
    # separated by colons
    def ctg_to_str(ctg_id, *args):
        return ':'.join([ctg_id, *map(lambda x: f'{x.get(ctg_id)}', args)])
    # returns string representing set of contigs
    def ctgs_to_str(ctg_ids, *args):
        return ','.join([ctg_to_str(ctg_id, *args) for ctg_id in ctg_ids])

    gc_bin, gc_term, gd_term, rd_term, path_contigs = 0, 0, 0, 0, []
    was_solved = model.getAttr('Status') == GRB.OPTIMAL

    gc_bins = pbm_input.get_gc_bins()
    bin1_gc_bin = gc_bins[bin1]
    bin1_gc = model.getVarByName('ctg_GC[s,{}]'.format(bin1_gc_bin)).Obj
    bin2_gc_bin = gc_bins[bin2]
    bin2_gc = model.getVarByName('ctg_GC[t,{}]'.format(bin2_gc_bin)).Obj

    pls_bins = pbm_input.get_pls_bins()
    bin1_ctgs = pls_bins.get_pls_content(with_mult=False)[bin1]
    bin2_ctgs = pls_bins.get_pls_content(with_mult=False)[bin2]
    overlap = set(bin1_ctgs) & set(bin2_ctgs) != set()
    ctg_lens = pbm_input.get_assembly().get_ctg_lens()
    bin1_ctgs_str = ctgs_to_str(bin1_ctgs, ctg_lens)
    bin2_ctgs_str = ctgs_to_str(bin2_ctgs, ctg_lens)
    flows = pls_bins.get_pls_copy_number()
    
    if was_solved:
        print('SUCCESS\tModel was solved')
        rd_term = model.getVarByName('d').X
        for v in model.getVars():
            if v.X > GUROBI_ZERO:
                if v.VarName[:2] == 'GC':
                    gc_bin = int(v.VarName[3])
                elif v.VarName[:3] == 'ctg':
                    gc_term += v.Obj
                elif v.VarName[0] == 'z':
                    gd_term += v.Obj
                    if v.VarName[2:-1] not in ['s', 't']:
                        path_contigs.append(v.VarName[2:-1])
    else:
        print('WARNING\tModel was not solved')                

    path_contigs = ctgs_to_str(path_contigs, ctg_lens)
    col_vals = [
        bin1, bin1_ctgs_str, bin1_gc_bin, bin1_gc, model.getVarByName('z[s]').Obj, flows[bin1],
        bin2, bin2_ctgs_str, bin2_gc_bin, bin2_gc, model.getVarByName('z[t]').Obj, flows[bin2],
        gc_bin, gc_term, gd_term, rd_term, path_contigs, int(not was_solved), overlap
    ]
    return col_vals

# Main merging function
def merging_all_pairs(
        sample,
        assembly_file,
        pls_scores_file,
        gc_intervals_file,
        pls_bins_file,
        source,
        model_sol_dir,
        out_tsv_file,
        threshold
):

    # Merges given bins and outputs list containing results
    def _merge_pair(pbm_input, bin1, bin2, sol_file, threshold):

        model = merging_model(pbm_input, bin1, bin2, threshold)
        model.setParam('OutputFlag', False)
        model.optimize()
        if model.getAttr('Status') == GRB.OPTIMAL:
            model.write(sol_file)
        line_data = merging_results(model, pbm_input, bin1, bin2)

        return [sample, *line_data]

    pbm_input = PBM_input(assembly_file, pls_scores_file, gc_intervals_file, pls_bins_file, source, gzipped=True)

    with open(out_tsv_file, 'w') as file:
        file.write('\t'.join(OUT_COLUMNS))
        pls_ids = pbm_input.get_pls_bins().get_pls_ids()
        if len(pls_ids) <= 1:
            print(f'{sample}.{source}: Less than 2 bins: no merging to consider')
        for bin1, bin2 in combinations(pls_ids, 2):
            sol_file = os.path.join(model_sol_dir, f'{sample}.{source}.{bin1}.{bin2}.sol')
            print(f'{sample}.{source}\t{bin1}.{bin2}\t{sol_file}')
            line = _merge_pair(pbm_input, bin1, bin2, sol_file, threshold)
            file.write('\n' + '\t'.join(['{}'.format(data) for data in line]))

# merges plasmid bins and removes contig repeats
def _flatten_bins(pred_bins, merged_ids):
    ctgs = {}
    for pred in merged_ids:
        ctgs.update(pred_bins[pred])
    return ','.join([f'{ctg[0]}:{ctg[1]}' for ctg in ctgs.items()])

# returns string to be written to post-merger plasmid bin TSV
def _flattened_bin_strs(pred_bins, merged_bins, flows):
    merged_ids, merged_ctgs = [], []
    for merger in merged_bins:
        merged_ids.append(str([pls for pls in merger]))
        merged_ctgs.append(_flatten_bins(pred_bins, merger))
    if flows:
        lines = zip(merged_ids, merged_ctgs, flows)
        return ['\t'.join([id, ctgs, str(rd)]) for id, ctgs, rd in lines]
    else:
        lines = zip(merged_ids, merged_ctgs)
        return ['\t'.join([id, ctgs]) for id, ctgs in lines]

# score based on how much read depth was lost merging the pair
def _combined_pbf_obj(results):
    post_obj = results[['MERGE_GC', 'MERGE_GD', 'MERGE_RD']].sum(axis=1)
    combined_obj = results[['BIN1_GC', 'BIN2_GC', 'BIN1_GD', 'BIN2_GD']].sum(axis=1) \
                                      + (results[['BIN1_RD', 'BIN2_RD']].sum(axis=1) / 2)
    return combined_obj - post_obj

# outputs merged plasmid bin file based on pair scoring
def merge_sample(
        assembly_file,
        pls_scores_file,
        gc_intervals_file,
        pls_bins_file,
        source,
        scored_sample_tsv,
        out_merger_file,
        score_func=_combined_pbf_obj,
        score_threshold=0
):
    
    results = pd.read_table(scored_sample_tsv)
    results['EDGE_WEIGHT'] = score_func(results)
    pbm_input = PBM_input(assembly_file, pls_scores_file, gc_intervals_file, pls_bins_file, source, gzipped=True)
    
    # constructing a graph where nodes are predicted bins and
    # edges are weighted by the pair's PlasMerge objective value
    G = nx.Graph()
    pls_bins = pbm_input.get_pls_bins()
    pls_content = pls_bins.get_pls_content(with_mult=True)
    G.add_nodes_from(pls_bins.get_pls_ids())
    nx.set_node_attributes(G, pls_bins.get_pls_copy_number(), name='rd')
    for index, row in results.iterrows():
        # edge only appears if PlasMerge found a solution
        if row['MILP_INFEASIBLE'] == 0:
            G.add_edge(row['BIN1'], row['BIN2'], weight=row['EDGE_WEIGHT'], rd=row['MERGE_RD'])
    # generate sample mergers by removing edges below weight threshold
    # and merging the resulting connected components
    above_thresh = [edge for edge in G.edges(data=True) if G[edge[0]][edge[1]]['weight'] >= score_threshold]
    H = nx.Graph()
    H.add_nodes_from(G.nodes(data=True))
    H.add_edges_from(above_thresh)
    merged_bins = [list(comp) for comp in nx.connected_components(H)]
    # if method is Plasbin-flow, take flow value to be min among merged bins
    flows = None
    if source == 'pbf':
        flows = [min(nx.get_edge_attributes(H.subgraph(comp), 'rd').values(), default=list(H.subgraph(comp).nodes(data='rd'))[0][1]) for comp in merged_bins]
    lines = _flattened_bin_strs(pls_content, merged_bins, flows)

    with open(out_merger_file, 'w') as file:
        if source == 'pbf':
            file.write('plasmid\tcontigs\tcopy_number')
        else:
            file.write('plasmid\tcontigs')
        file.writelines(['\n' + line for line in lines])

if __name__ == "__main__":
    import os
   
    samples = ['SAMN32247302', 'SAMN32247345', 'SAMN32247425', 'SAMN32247519', 'SAMN32247522']
    #samples = ['SAMN32247519']
    root = os.path.normpath('../test')
    gc_intervals_file = os.path.join(root, 'gc_intervals.txt')
    sources = ['gt', 'mob', 'gp', 'pbf']
    threshold = 0.05
   
    for sample in samples:
        print(f'SAMPLE: {sample}')
       
        gfa_file = os.path.join(root, 'gfas', f'{sample}.assembly.gfa.gz')
        pls_scores_file = os.path.join(root, 'scores', f'{sample}.scores.tsv')
        model_sol_dir = os.path.join(root, 'results', 'model')

        for source in sources:
            print(f'\tSOURCE: {source}')
           
            pls_bins_file = os.path.join(root, 'pls_bins', f'{sample}.{source}.tsv')
            out_tsv_file = os.path.join(root, 'results', source, f'{sample}.{source}.tsv')
            out_merger_file = os.path.join(root, 'results', source, f'{sample}.{source}.mergers.txt')
           
            merging_all_pairs(
                sample,
                gfa_file,
                pls_scores_file,
                gc_intervals_file,
                pls_bins_file,
                source,
                model_sol_dir,
                out_tsv_file,
                threshold
            )

            merge_sample(
                gfa_file,
                pls_scores_file,
                gc_intervals_file,
                pls_bins_file,
                source,
                out_tsv_file,
                out_merger_file
            )
