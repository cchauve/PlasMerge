#!/usr/bin/python

from __future__ import division
import argparse
import os
import pandas as pd
import io_utils

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

def eval_bins(sid, sample_gt, sample_pbf, sample_assembly, recall, precision, len_th):
	gt_plasmids = sample_gt.plasmids
	pbf_bins = sample_pbf.plasmids
	ctg_details = sample_assembly.ctg_dict

	for pls_id in gt_plasmids:
		ref_pls = sid+'_GT_'+pls_id
		recall[ref_pls] = {'wtd': {}, 'unwtd': {}}
		for eval_type in ['wtd', 'unwtd']: 
			recall[ref_pls][eval_type] = {'Val': 0, 'Bin': 'NA', 'Common': 0, 'Total': 0}		
		nref_ctgs = 0
		lenref_ctgs = 0
		for ctg in gt_plasmids[pls_id]:
			ctg_len = ctg_details[ctg[0]].length
			if ctg_len >= len_th:
				nref_ctgs += 1
				lenref_ctgs += ctg_len
		recall[ref_pls]['unwtd']['Total'] = nref_ctgs
		recall[ref_pls]['wtd']['Total'] = lenref_ctgs					
		for bin_id in pbf_bins.plasmids:
			pred_pls = sid + '_' + bin_id
			common_ctgs = set([ctg[0] for ctg in pbf_bins.plasmids[bin_id]]).intersection(set([ctg[0] for ctg in gt_plasmids[pls_id]]))
			ncommon_ctgs = 0
			lencommon_ctgs = 0
			for ctg in common_ctgs:
				ctg_len = ctg_details[ctg].length
				if ctg_len >= len_th:
					ncommon_ctgs += 1
					lencommon_ctgs += ctg_len
			n_rec, len_rec = 0, 0
			if nref_ctgs >= 1:
				n_rec = ncommon_ctgs / nref_ctgs
				len_rec = lencommon_ctgs / lenref_ctgs	
			if n_rec > recall[ref_pls]['unwtd']['Val']:
				recall[ref_pls]['unwtd'] = {'Val': n_rec, 'Bin': pred_pls, 'Common': ncommon_ctgs, 'Total': nref_ctgs}	
			if len_rec > recall[ref_pls]['wtd']['Val']:
				recall[ref_pls]['wtd'] = {'Val': len_rec, 'Bin': pred_pls, 'Common': lencommon_ctgs, 'Total': lenref_ctgs}	

	for bin_id in pbf_bins.plasmids:
		pred_pls = sid + '_' + bin_id
		precision[pred_pls] = {'wtd': {}, 'unwtd': {}}
		for eval_type in ['wtd', 'unwtd']: 
			precision[pred_pls][eval_type] = {'Val': 0, 'Ref': 'NA', 'Common': 0, 'Total': 0}
		npred_ctgs = 0
		lenpred_ctgs = 0	
		for ctg in pbf_bins.plasmids[bin_id]:
			ctg_len = ctg_details[ctg[0]].length
			if ctg_len >= len_th:
				npred_ctgs += 1
				lenpred_ctgs += ctg_len
		precision[pred_pls]['unwtd']['Total'] = npred_ctgs
		precision[pred_pls]['wtd']['Total'] = lenpred_ctgs
		for pls_id in gt_plasmids:
			ref_pls = sid+'_GT_'+pls_id
			common_ctgs = set([ctg[0] for ctg in gt_plasmids[pls_id]]).intersection(set([ctg[0] for ctg in pbf_bins.plasmids[bin_id]]))
			ncommon_ctgs = 0
			lencommon_ctgs = 0
			for ctg in common_ctgs:
				ctg_len = ctg_details[ctg].length
				if ctg_len >= len_th:
					ncommon_ctgs += 1
					lencommon_ctgs += ctg_len
			n_prec, len_prec = 0, 0
			if npred_ctgs >= 1:
				n_prec = ncommon_ctgs / npred_ctgs
				len_prec = lencommon_ctgs / lenpred_ctgs
			if n_prec > precision[pred_pls]['unwtd']['Val']:
				precision[pred_pls]['unwtd'] = {'Val': n_prec, 'Ref': ref_pls, 'Common': ncommon_ctgs, 'Total': npred_ctgs}
			if len_prec > precision[pred_pls]['wtd']['Val']:
				precision[pred_pls]['wtd'] = {'Val': len_prec, 'Ref': ref_pls, 'Common': lencommon_ctgs, 'Total': lenpred_ctgs}

	return recall, precision
#---------------------------------------------------------

#---------------------------------------------------------
#Listing all the samples
EFAECIUM = [15,16]
ECOLI = [18,19,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,52]
KPNEUMONIAE = [62,63,64,65,66,76,85,86,87]
OTHER = [1,5,55,56,102,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,129,133]
TEST_SAMPLES = EFAECIUM + ECOLI + KPNEUMONIAE + OTHER

def calculate_pbf_stats(merge_type=None, merge_thr=None, aggregate=True):
    
	pls_all = {}
	pbf_all = {}
	assembly = {}
	recall = {}
	precision = {}
    
	for sample in TEST_SAMPLES:
		sid = str(sample)
		ground_truth = '../../../organized_data/ground_truth/sample_'+sid+'/ctg_pls.tsv'	#PATH TO GROUND TRUTH
		gc_int_file = '../../../organized_data/gc_intervals.txt'	#PATH TO GC INTERVALS FILE
		assembly_file = '../../../plasbin_flow_input/sample_'+sid+'/assembly.gfa'	#PATH TO ASSEMBLY FILE
		if merge_thr != None:
			pbf_file = '../../../organized_data/post_processing_mergers_'+merge_type+'/sample_'+sid+'/'+str(merge_thr)+'/plasbin_flow_bins.out'	#PATH TO PLASBIN-FLOW BINS
		else:
			pbf_file = '../../../organized_data/output/sample_'+sid+'/plasbin_flow_bins.out'

		#For the sample:
		#Computing recall for reference plasmids and
		#Computing precision for predicted bins	
		#Updating recall and precision dictionaries
		if os.path.isfile(ground_truth):
			assembly[sid] = io_utils.Assembly(assembly_file)
			pls_all[sid] = io_utils.Plasmids('ground_truth', ground_truth, 0.95)
			pbf_all[sid] = io_utils.PBF_output('plasbin_flow', pbf_file, gc_int_file)

			eval_bins(sid, pls_all[sid], pbf_all[sid], assembly[sid], recall, precision, 0)

	print('Stats for threshold', merge_thr)

	if aggregate:
		#Computing average precision (unweighted and weighted)
		prec_n_common, prec_len_common = 0, 0
		prec_n_total, prec_len_total = 0, 0
		for p in precision:
			prec_n_common += precision[p]['unwtd']['Common']
			prec_n_total += precision[p]['unwtd']['Total']
			prec_len_common += precision[p]['wtd']['Common']
			prec_len_total += precision[p]['wtd']['Total']	
		ovr_n_prec = prec_n_common / prec_n_total
		ovr_len_prec = prec_len_common / prec_len_total
		print('Precision (unweighted, weighted):', ovr_n_prec, ovr_len_prec)

		#Computing average recall (unweighted and weighted)
		rec_n_common, rec_len_common = 0, 0
		rec_n_total, rec_len_total = 0, 0
		for p in recall:
			rec_n_common += recall[p]['unwtd']['Common']
			rec_n_total += recall[p]['unwtd']['Total']
			rec_len_common += recall[p]['wtd']['Common']
			rec_len_total += recall[p]['wtd']['Total']	
		ovr_n_rec = rec_n_common / rec_n_total
		ovr_len_rec = rec_len_common / rec_len_total
		print('Recall (unweighted, weighted):', ovr_n_rec, ovr_len_rec)

		#Computing F1 score (unweighted and weighted)
		ovr_n_f1 = 2*ovr_n_prec*ovr_n_rec / (ovr_n_prec + ovr_n_rec)
		ovr_len_f1 = 2*ovr_len_prec*ovr_len_rec / (ovr_len_prec + ovr_len_rec)
		print('F1 (unweighted, weighted):', ovr_n_f1, ovr_len_f1)
		return ovr_n_prec, ovr_len_prec, ovr_n_rec, ovr_len_rec, ovr_n_f1, ovr_len_f1
	else:
		prec_n, prec_len = {}, {}
		for p in precision:
			prec_n[p] = precision[p]['unwtd']['Common'] / precision[p]['unwtd']['Total']
			prec_len[p] = precision[p]['wtd']['Common'] / precision[p]['wtd']['Total']
		rec_n, rec_len = {}, {}
		for p in recall:
			rec_n = recall[p]['unwtd']['Common'] / recall[p]['unwtd']['Total']
			rec_len[p] = recall[p]['wtd']['Common'] / recall[p]['wtd']['Total']
		return prec_n, prec_len, rec_n, rec_len