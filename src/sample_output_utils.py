import io_utils as io
import networkx as nx
import numpy as np
import os

# compiles assembly, predicted bins and ground truth for each sample
# (requires that all sample folders are in the same directory 'sample_path'
# and similarly for ground truth)
def compile_pbf_results(gt_path, pred_path, gc_int_path):
    
    assembly_dict = {}
    pred_bin_dict = {}
    gt_dict = {}
    
    for sample in os.listdir(pred_path):
        assembly_file = gt_path + '/' + sample + '/assembly.gfa'
        gt_map = gt_path + '/' + sample + '/ctg_pls.tsv'
        bins_file = pred_path + '/' + sample + '/plasbin_flow_bins.out'
        assembly_dict[sample] = io.Assembly(assembly_file)
        gt_dict[sample] = io.Plasmids('ground_truth', gt_map)
        pred_bin_dict[sample] = io.PBF_output('plasbin_flow', bins_file, gc_int_path)
    
    return assembly_dict, gt_dict, pred_bin_dict 

# outputs formatted string representing a plasmid from its contig IDs
def __plasmid_str(ctg_ids, ctg_dict):
    
    out_str = '('
    
    for ctg in ctg_ids:
        out_str += ctg + ':' + str(ctg_dict[ctg].get_len())
        if ctg != ctg_ids[-1]:
            out_str += ','
    
    out_str += ')'
    
    return out_str
    
# generates string for which predicted bins are properly contained in
# ground truth plasmids for a given sample
def __sample_pred_containment(sample, assembly, pred_bins, gt):
    
    out = ''
    
    for pls in list(gt.get_pls_ids()):
        contained_pred = []
        contained_pred_ids = []
        out += sample + '\t' + pls + '\t'
        ctg_ids = [ctg[0] for ctg in gt.plasmids[pls]]
        out += __plasmid_str(ctg_ids, assembly.ctg_dict) + '\t'
        
        for pred_pls in list(pred_bins.plasmids.get_pls_ids()):
            pred_ids = [ctg[0] for ctg in pred_bins.plasmids.plasmids[pred_pls]]
            if (set(pred_ids) & set(ctg_ids)) == set(pred_ids):
                contained_pred_ids.append(pred_pls)
                contained_pred.append(__plasmid_str(pred_ids, assembly.ctg_dict))
        
        if contained_pred:
            out += '('
            for pls_id in contained_pred_ids:
                out += pls_id
                if pls_id != contained_pred_ids[-1]:
                    out += ','
            out += ')\t'
            out += '('
            for pls_str in contained_pred:
                out += pls_str
                if pls_str != contained_pred[-1]:
                    out += ';'
            out += ')'
                    
        out += '\n'
    
    return out

# as above but for a dictionary of samples; optionally writes output file
def pred_containment(assembly_dict, gt_dict, bin_dict, outpath=None):
    
    out = 'SAMPLE\tPLASMID\tGT\tPRED_IDS\tPRED_BINS\n'
    for sample in list(assembly_dict.keys()):
        out += __sample_pred_containment(sample, assembly_dict[sample], \
                                               bin_dict[sample], gt_dict[sample])
    if outpath:
        with open(outpath + '/splitting.tsv', 'w') as tsv_file:
            tsv_file.write(out)
    else:
        return out

# plots subgraph induced by nodes from predicted plasmids (optionally including
# other nodes within the given radius), with nodes
# coloured based on predicted bin containment
# (requires that nodes of G have 'pred_bins' and 'gt_bin' attributes, which
# are lists of predicted and ground truth bins the node is in)
def view_plasmids_subgraph(G, radius=0):
    
    plasmids = set()
    plasmid_nodes = []
    node_labels = {}
    node_colours = []
    node_label_conv = set()
    extended_node_set = set()
    extra_nodes = set()
    
    for node in G.nodes():
        if G.nodes[node]['pred_bins']:
            plasmids.add(G.nodes[node]['pred_bins'][0])    
            plasmid_nodes.append(node)
            node_labels[node] = G.nodes[node]['gt_bin']
            node_label_conv.update(G.nodes[node]['gt_bin'])
    
    node_label_conv = list(node_label_conv)
    plasmids = list(plasmids) 
    node_label_conv = {plasmid: i for plasmid, i in zip(node_label_conv, range(20))}
    colours = {plasmid: i for plasmid, i in zip(plasmids, np.arange(0.30, 1, 0.03))}
    print(colours)

    if radius > 0:
        for node in plasmid_nodes:
            within_radius = nx.single_source_shortest_path_length(G, node, radius).keys()
            extended_node_set = extended_node_set.union(within_radius)
        extra_nodes = extended_node_set.difference(plasmid_nodes)
        H = nx.induced_subgraph(G, extended_node_set)
    else:
        H = nx.induced_subgraph(G, plasmid_nodes)
    
    for node in H.nodes():
        if node in extra_nodes:
            node_labels[node] = ""
            node_colours.append(1.0)
        else:
            node_labels[node] = [node_label_conv[plasmid] for plasmid in node_labels[node]]
            node_colours.append(colours[H.nodes[node]['pred_bins'][0]]) 
    
    try:
        pos = nx.planar_layout(H, scale=2)
    except:
        pos = nx.spring_layout(H, k=0.6, scale=2)
        
    nx.draw_networkx(H, pos, labels=node_labels, node_size=150, node_color=node_colours, font_color='cyan', cmap='magma')
    