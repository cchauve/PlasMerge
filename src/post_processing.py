import os
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import io_utils as io

# sums up GC and gene density terms for contigs in bin_ctgs
# into single node with name node
def collapse_obj_terms(gd, gc_probs, gc_ints, bin_ctgs, node):
    
    merged_gc = {}
    merged_gd = {}
    
    for ctg in bin_ctgs:
        max_prob = max(gc_probs[ctg])
        if merged_gd.get(node) != None:
            for i in gc_ints:
                merged_gc[(node, i)] += gc_probs[ctg][i] - max_prob
            merged_gd[node] += gd[ctg] - 0.5
        else:
            for i in gc_ints:
                merged_gc[(node, i)] = gc_probs[ctg][i] - max_prob
            merged_gd[node] = gd[ctg] - 0.5

    return merged_gc, merged_gd
        
# computes and returns objective function terms for each node 
# in modified graph (after contracting both bins into one node)
# to be used in post-processing MILP
def obj_terms(gd, gc_probs, gc_ints, b1_ctgs, b2_ctgs):
    
    b1_gc, b1_gd = collapse_obj_terms(gd, gc_probs, gc_ints, b1_ctgs, 's')
    b2_gc, b2_gd = collapse_obj_terms(gd, gc_probs, gc_ints, b2_ctgs, 't')
    b1_gc.update(b2_gc)
    b1_gd.update(b2_gd)
    
    for ctg in set(gd.keys()).difference(set(b1_ctgs) | set(b2_ctgs)):
        max_prob = max(gc_probs[ctg])
        for i in gc_ints:
            b1_gc[(ctg, i)] = gc_probs[ctg][i] - max_prob
        b1_gd[ctg] = gd[ctg] - 0.5
      
    return b1_gc, b1_gd

# Removes read_depth used up by other bins
def remove_plasmid_rd(rd, flows, pred_bins, pred_ids, bin1, bin2):
    
    for key in pred_ids:
        if key in [bin1, bin2]:
            continue
        bin_mults = dict(pred_bins.get_pls(key))
        for ctg in bin_mults:
            if rd.get(ctg, None) == None:
                continue
            else:
                rd[ctg] = rd[ctg] - (bin_mults[ctg] * flows[key])
                
# Removing vertices that now have read-depth less than threshold                
def remove_low_rd(gc_pen, gd_pen, rd, G, gc_ints, threshold):
    
    to_remove = [node for node in rd if rd[node] < threshold]
    #print('to remove:', to_remove)
    for node in to_remove:
        for i in gc_ints:
            gc_pen.pop((node, i))
        gd_pen.pop(node)
        rd.pop(node)
        G.remove_node(node)

    
# contracting bins down to single vertices in networkx graph
def contracted_bins(graph, b1_ctgs, b2_ctgs):
    
    H = nx.DiGraph(graph)
    b1_set = set(b1_ctgs)
    b2_set = set(b2_ctgs)
    shared_ctgs = list(b1_set & b2_set)
    b1_diff = list(b1_set - b2_set)
    b2_diff = list(b2_set - b1_set)
    for ctg in b1_diff[1:]:
        nx.contracted_nodes(H, b1_diff[0], ctg, self_loops=False, copy=False)
    for ctg in b2_diff[1:]:
        nx.contracted_nodes(H, b2_diff[0], ctg, self_loops=False, copy=False)
    H = nx.relabel_nodes(H, {b1_diff[0]:'s', b2_diff[0]:'t'})
    if shared_ctgs:
        for ctg in shared_ctgs[1:]:
            nx.contracted_nodes(H, shared_ctgs[0], ctg, self_loops=False, copy=False)
        shared_edges = nx.edge_boundary(H, [shared_ctgs[0]], set(H.nodes()) - {'s', 't'})
        for source, nbr in shared_edges:
            H.add_edges_from([('s', nbr), (nbr, 's'), ('t', nbr), (nbr, 't')])
        H.add_edges_from([('s', 't'), ('t', 's')])
        H.remove_node(shared_ctgs[0])

    return H

# Creates MILP gurobipy object to be optimized
def generate_model(H, gc_pen, gd_pen, rd, gc_ints):
    
    model = gp.Model("model")
    model.setAttr('ModelSense', GRB.MAXIMIZE)
    edge_list = list(H.edges())
    node_list = list(H.nodes())

    x = model.addVars(edge_list, vtype=GRB.BINARY, name="x")
    q = model.addVars(edge_list, lb=0, vtype=GRB.CONTINUOUS, name="q")
    z = model.addVars(node_list, obj=gd_pen, vtype=GRB.BINARY, name="z")
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
def post_processing_model(sample, bin1, bin2, input_dir, output_dir, gc_file, threshold):
    
    pbf_input = io.PBF_input(os.path.join(input_dir, sample, 'assembly.gfa'), gc_file, \
                             os.path.join(input_dir, sample, 'filtered_genes_to_contigs.csv'), 0.58, 2650)
    pbf_output = io.PBF_output('plasbin_flow', os.path.join(output_dir, sample, 'plasbin_flow_bins.out'), gc_file)
    pred_bins = pbf_output.get_pls()
    pred_ids = pred_bins.get_pls_ids()
    flows = pbf_output.flow
    b1_ctgs = [ctg[0] for ctg in pred_bins.get_pls(bin1)]
    b2_ctgs = [ctg[0] for ctg in pred_bins.get_pls(bin2)]
    G = pbf_input.to_graph().to_directed()
    gd = pbf_input.get_gd()
    gc_probs = pbf_input.get_gc_probs()
    rd = {}
    gc_ints = range(len(pbf_input.get_gc_intervals()) - 1)

    gc_pen, gd_pen = obj_terms(gd, gc_probs, gc_ints, b1_ctgs, b2_ctgs)
    for ctg in set(gd.keys()).difference(set(b1_ctgs) | set(b2_ctgs)):
            rd[ctg] = pbf_input.get_ctg(ctg).get_rd()
    rd['s'] = float(flows[bin1])
    rd['t'] = float(flows[bin2])
    remove_plasmid_rd(rd, flows, pred_bins, pred_ids, bin1, bin2)
    H = contracted_bins(G, b1_ctgs, b2_ctgs)
    if threshold == 'min':
        remove_low_rd(gc_pen, gd_pen, rd, H, gc_ints, min(rd['s'], rd['t']))
    else:
        remove_low_rd(gc_pen, gd_pen, rd, H, gc_ints, threshold)
    
    return generate_model(H, gc_pen, gd_pen, rd, gc_ints)

# returns whether lst1 is contained in lst2
def list_contained(lst1, lst2):
    return set(lst1) & set(lst2) == set(lst1)


# returns predicted bins with those properly contained
# in another removed
def filtered_bins(sample, output_file, gc_file):

    pbf_out = io.PBF_output('plasbin_flow', output_file, gc_file)
    pred_bins = pbf_out.get_pred_bins(with_mults=False)
    contained_mask = {}
    for pred_id in pred_bins.keys():
        contained_mask[pred_id] = any(list_contained(pred_bins[pred_id], other) \
                                      for other in pred_bins.keys() if other != pred_id)
    
    return [key for key, val in pred_bins.items() if val]

########################### DISPLAYING/IO ###########################

# prints some solution details from solution file
def print_post_processing_sol_file(sol_file):
    
    nz_flows = {}
    rd = 0
    with open(sol_file, 'r') as file:
        line = next(file, None)
        line = next(file, None)
        obj = float(line[:-1].split(sep=' ')[-1])
        line = next(file, None)
        while line:
            lst = line[:-1].split(sep=' ')
            if float(lst[-1]) < 0.01:
                line = next(file, None)
                continue
            if lst[0][0] == 'q':
                nz_flows[lst[0]] = float(lst[-1])
            elif lst[0][0] == 'd':
                rd = float(lst[-1])
            line = next(file, None)

    print('Optimal path (as defined by flow values):')
    for flow in sorted(nz_flows, key=nz_flows.get, reverse=True):
        print(flow, nz_flows[flow])

    print('\nMinimum read depth:', rd)
    print('Objective value:', obj)