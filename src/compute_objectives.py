#!/usr/bin/python

from __future__ import division
import numpy as np
#import argparse
import io_utils

#Functions to compute objective values
def compute_best_GC_bin(ctg_list, pbf_input):
	'''
	Computes the GC bin with the highest frequency 
	of having highest GC prob for a given list of contigs, 
	'''
	best_bin = []
	for ctg in ctg_list:
		ctg_id = ctg[0]
		best_bin.append(np.argmax(pbf_input.GC_probs[ctg_id]))
	return max(set(best_bin), key = best_bin.count)	+ 1

def compute_GC_pens(ctg_list, pbf_input, GC_bin):
	'''
	Computes the GC penalties per contig and the total GC penalty 
	associated with a given list of contigs and GC_bin
	'''
	GC_pens, total_GC_pen = {}, 0
	#GC_pens = []
	for ctg in ctg_list:
		ctg_id = ctg[0]
		max_prob = max(pbf_input.GC_probs[ctg_id])
		GC_pens[ctg_id] = [p - max_prob for p in pbf_input.GC_probs[ctg_id]][GC_bin - 1]
		total_GC_pen += GC_pens[ctg_id]
	return GC_pens, total_GC_pen

def compute_total_gd(ctg_list, pbf_input, offset = 0.5):
	'''
	Computes the total gene density for a given list of contig
	'''
	total_gd = 0
	for ctg in ctg_list:
		ctg_id = ctg[0]
		total_gd += pbf_input.gd[ctg_id] - offset
	return total_gd
