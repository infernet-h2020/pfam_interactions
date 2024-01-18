import os
import sys
import gzip
import time
import codecs
import datetime
import subprocess
import multiprocessing
import argparse
import pickle
import Bio.PDB
from sklearn.cluster import DBSCAN
from sklearn import metrics
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from Bio.SubsMat import MatrixInfo as matlist
import numpy as np
from collections import Mapping, Container
from sys import getsizeof
import shutil
import re
#import objgraph


def linear_response(d, a=8, b=8):
	if b<a:
		print("Error: must be a <= b")
		exit(1)
	if d<=a:
		return 1
	elif d>=b:
		return 0
	else:
		return 1 - (d-a)/(b-a)
	

def from3to1(resname, include_MSEs=False):
	f3t1 = {'ALA' : 'A',
	        'ARG' : 'R',
	        'ASN' : 'N',
	        'ASP' : 'D',
	        'CYS' : 'C',
	        'GLN' : 'Q',
	        'GLU' : 'E',
	        'GLY' : 'G',
	        'HIS' : 'H',
	        'ILE' : 'I',
	        'LEU' : 'L',
	        'LYS' : 'K',
	        'MET' : 'M',
	        'PHE' : 'F',
	        'PRO' : 'P',
	        'SER' : 'S',
	        'THR' : 'T',
	        'TRP' : 'W',
	        'TYR' : 'Y',
	        'VAL' : 'V'}

	if resname in list(f3t1.keys()):
		return f3t1[resname]
	elif include_MSEs and resname == 'MSE':
		return 'M'
	else:
		return 'X'

def query_yes_no(question, default="yes"):
	"""Ask a yes/no question via input() and return their answer.

	"question" is a string that is presented to the user.
	"default" is the presumed answer if the user just hits <Enter>.
	It must be "yes" (the default), "no" or None (meaning
	an answer is required of the user).

	The "answer" return value is True for "yes" or False for "no".
	"""
	valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
	if default is None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '%s'" % default)

	while True:
		sys.stdout.write(question + prompt)
		choice = input().lower()
		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")


def string_is_float(s):
	if not type(s) == str:
		print("ERROR: input type {0}, while str expected".format(type(s)))
		exit(1)
	try:
		float(s)
		return True
	except ValueError:
		return False


def string_is_int(s):
	if not type(s) == str:
		print("ERROR: input type {0}, while str expected".format(type(s)))
		exit(1)
	try:
		int(s)
		return True
	except ValueError:
		return False


def string_isnot_dict(element):
	inquotes = ""
	for c in element:
		if c == '"' or c == "'":
			if c == inquotes:
				inquotes = ""
			elif not inquotes:
				inquotes = c
		elif c == ":" and not inquotes:
			return False
	return True


def string_is_properstring(element):
	if len(element) < 2 or (not ((element[0] == "'" and element[-1] == "'" ) or (element[0] == '"' and element[-1] == '"'))):
		return False
	for nc, c in enumerate(element[1:-1]):
		if c == element[0] and (nc == 0 or element[nc-1] != "\\" ):
			return False
	return True


def string_decode_element(element, simple_version=False, is_immutable=False, permissive=False):
	"""WARNING: Cannot decode following types:
	decimal, complex, range, bytes, bytearrary,
	and any mutable or immutable user-defined class
	or non-standard library class (e.g. np.array())
	"""
	global log_filename
	this_name = string_decode_element.__name__

	nop_err = ()
	element = element.strip()
	if simple_version:
		if element == '':
			return ''
		elif element == "None":
			return None
		elif element == "False":
			return False
		elif element == "True":
			return True
		elif string_is_int(element):
			return int(element)
		elif string_is_float(element):
			return float(element)
		else:
			return element

	if string_is_properstring(element):
		if len(element) == 2:
			return ""
		else:
			return element[1:-1]
	elif element == "None":
		return None
	elif element == "False":
		return False
	elif element == "True":
		return True
	elif string_is_int(element):
		return int(element)
	elif string_is_float(element):
		return float(element)
	elif (len(element) > 1 and element[0] == "(" and element[-1] == ")") or (len(element) > 6 and element[:6] == "tuple(" and element[-1] == ")"):
		return string_decode_list(element, is_tuple=True)
	elif len(element) > 10 and element[:10] == "frozenset(" and element[-1] == ")":
		return string_decode_list(element, is_fset=True)
	elif not is_immutable:
		if (len(element) > 1 and element[0] == "{" and element[-1] == "}" and string_isnot_dict(element)) or (len(element) > 4 and element[:4] == "set(" and element[-1] == ")"):
			return string_decode_list(element, is_set=True)
		elif (len(element) > 1 and element[0] == "[" and element[-1] == "]") or (len(element) > 5 and element[:5] == "list(" and element[-1] == ")"):
			return string_decode_list(element)
		elif (len(element) > 1 and element[0] == "{" and element[-1] == "}") or (len(element) > 5 and element[:5] == "dict(" and element[-1] == ")"):
			return string_decode_dict(element)
		elif permissive:
			return element
		else:
			nop_err = "ERROR: could not process element {0}".format(element)
	elif permissive:
		return element
	else:
		nop_err = "ERROR: could not process immutable element {0}".format(element)

	if nop_err:
		print(nop_err)
		exit(1)


def convindex_uniprot_dca__format_stockholm(line):
	"""Converts indices from uniprot to DCA for a protein sequence in Stockholm format

	# Input
	(string) line: sequence in Stockholm format
	Stockholm format: 2 fields (WARNING: mind the format of field 0)
	   field 0: {UniProt_accession}/{UniProt_initial_resID}-{UniProt_final_resID} 
	   field 1: sequence (uppercase: matches, lowercase: insertions, dash: deletions, dot: placeholder)

	# Output
	(list) conversion: list of lists of length 2
	   entry 0: header pair [(int), (int)] [UniProt_initial_resID, UniProt_final_resID]
	   entries 1 to N: conversion pairs [(int), (int)] [DCA_resID, UniProt_resID]
	(dict) uniprot_id2name: converts UniProt positions (resIDs) into the corresponding amino acid
	   keys: (int) UniProt resID
	   values: (char) UniProt residue types (1-letter)
	(int) dca_resid: length of the DCA model
	"""

	# Parse line
	fields = line.split()
	if len(fields) != 2:
		print("ERROR: line must contain 2 fields", line)
		exit(1)
	initial_uniprot_resid = int(fields[0].split('/')[1].split('-')[0])
	final_uniprot_resid = int(fields[0].split('/')[1].split('-')[1])
	seqaln = fields[1]

	# Compute conversion and UniProt chain
	conversion = [[initial_uniprot_resid, final_uniprot_resid]]
	uniprot_id2name = {}

	dca_resid = 1
	uniprot_resid = initial_uniprot_resid
	for c in seqaln:
		if c == '-':
			dca_resid += 1
		elif c == '.':
			pass
		elif c == c.lower():
			uniprot_id2name[uniprot_resid] = c.upper()
			uniprot_resid += 1
		elif c == c.upper():
			conversion.append([dca_resid, uniprot_resid])
			uniprot_id2name[uniprot_resid] = c
			dca_resid += 1
			uniprot_resid += 1


	return conversion, uniprot_id2name, dca_resid-1


def download_pfam_files(pfam_acc, folder, msa_type, version, only_name=False):
	folder = folder + pfam_acc[:4] + '/'
	return folder + "/{0}_{1}_v{2}.stockholm.gz".format(pfam_acc, msa_type, version)
	if version != 32:
		print(version)
		print("THE version OPTION IS NOT YET IMPLEMENTED")
		exit(1)
	else:	# Assumes you need the latest version
		if msa_type == 'hmm':
			if not only_name:
				subprocess.run(["wget", "https://pfam-legacy.xfam.org/family/{0}/hmm".format(pfam_acc), "-O", "{1}/{0}.hmm".format(pfam_acc, folder)], stdout=open("/dev/null", 'w'))
			return folder + "/{0}.hmm".format(pfam_acc)
		else:
#			print("wget https://pfam.xfam.org/family/{0}/alignment/{1} -O {2}/{0}_{1}_v{3}.stockholm".format(pfam_acc, msa_type, folder, version), only_name)
			if not only_name:
#				print("wget https://pfam.xfam.org/family/{0}/alignment/{1} -O {2}/{0}_{1}_v{3}.stockholm".format(pfam_acc, msa_type, folder, version))
				subprocess.run(["wget", "https://pfam-legacy.xfam.org/family/{0}/alignment/{1}".format(pfam_acc, msa_type), "-O", "{2}/{0}_{1}_v{3}.stockholm".format(pfam_acc, msa_type, folder, version)], stdout=open("/dev/null", 'w'))
				if os.path.exists("{2}/{0}_{1}_v{3}.stockholm".format(pfam_acc, msa_type, folder, version)):
					print("ERROR: Could not download the file https://pfam-legacy.xfam.org/family/{0}/alignment/{1}\nThe algorithm cannot proceed. Please download it and change its path in {2}/{0}_{1}.stockholm".format(pfam_acc, msa_type, folder))
					exit(1)
			return folder + "/{0}_{1}_v{2}.stockholm.gz".format(pfam_acc, msa_type, version)


def download_pdb(pdbname, output_folder):
	print("wget", "https://files.rcsb.org/download/{0}.pdb".format(pdbname.upper()), "-O", "{0}/{1}.pdb")
	subprocess.run(["wget", "https://files.rcsb.org/download/{0}.pdb".format(pdbname.upper()), "-O", "{0}/{1}.pdb".format(output_folder, pdbname.lower())], stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))


def calculate_asymm_seqid(alignment):
	ntot1 = 0
	ntot2 = 0
	naln = 0
	for na in range(len(alignment[0])):
		if alignment[0][na] != '-':
			ntot1 += 1
		if alignment[1][na] != '-':
			ntot2 += 1
		if alignment[0][na] == alignment[1][na]:
			naln += 1

	if ntot1 == 0:
		seqID1 = 0
	else:
		seqID1 = naln/ntot1
	if ntot2 == 0:
		seqID2 = 0
	else:
		seqID2 = naln/ntot2

	return seqID1, seqID2


def instrinsic_dimension_clustering(X):
	pass
	# Call intrinsic dimension

	# Call iterative k-means

	# Take the one that maximizes the silhouette score and return its labels


def dbscan(X):
	# #############################################################################
	# Compute DBSCAN
	db = DBSCAN(eps=0.1, min_samples=10).fit(X)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_
	
	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	n_noise_ = list(labels).count(-1)
	
	print('Estimated number of clusters: %d' % n_clusters_)
	print('Estimated number of noise points: %d' % n_noise_)
#	print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#	print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#	print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#	print("Adjusted Rand Index: %0.3f"
#	      % metrics.adjusted_rand_score(labels_true, labels))
#	print("Adjusted Mutual Information: %0.3f"
#	      % metrics.adjusted_mutual_info_score(labels_true, labels))
#	print("Silhouette Coefficient: %0.3f"
#	      % metrics.silhouette_score(X, labels))
	
	# #############################################################################
	# Plot result
	
	# Black removed and is used for noise instead.
	unique_labels = set(labels)
	colors = [plt.cm.Spectral(each)
	          for each in np.linspace(0, 1, len(unique_labels))]
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = [0, 0, 0, 1]
	
	    class_member_mask = (labels == k)
	
	    xy = X[class_member_mask & core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=14)
	
	    xy = X[class_member_mask & ~core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)
	
	plt.title('Estimated number of clusters: %d' % n_clusters_)
	plt.savefig("fig1.png")

	with open("m1.txt", 'w') as m_file:
		N = len(X[0])
		for i in range(N):
			for j in range(N):
				m_file.write("{0} {1} {2}\n".format(i, j, X[i,j]))


def deep_getsizeof(o, ids):
	"""Find the memory footprint of a Python object
	
	This is a recursive function that drills down a Python object graph
	like a dictionary holding nested dictionaries with lists of lists
	and tuples and sets.
	
	The sys.getsizeof function does a shallow size of only. It counts each
	object inside a container as pointer only regardless of how big it
	really is.

	:param o: the object
	:param ids:
	:return:

	Call: deep_getsizeof(object_to_be_evaluated, set())
	"""
	d = deep_getsizeof
	if id(o) in ids:
		return 0

	r = getsizeof(o)
	ids.add(id(o))

	if isinstance(o, str):
		return r

	if isinstance(o, Mapping):
		return r + sum(d(k, ids) + d(v, ids) for k, v in o.items())

	if isinstance(o, Container):
		return r + sum(d(x, ids) for x in o)
	
	return r


def check_pdb_quality(pdb_path, resolution_threshold):
	with open(pdb_path) as f:
		for line in f:
			if line.startswith("EXPDTA"):
				if line.split()[1] != "X-RAY":
					print(pdb_path, "not an X-ray structure")
					return True#False
			if line.startswith("REMARK   2 RESOLUTION"):
				if float(line.split()[3]) > resolution_threshold:
					print(pdb_path, "not enough resolution:", line.split()[3], ">", resolution_threshold)
					return False
				print("PDB good!")
				return True

def extract_sequence_from_pdb(pdb_path, uniprot_spec=False):
	pdb_seqs = {}
	uniprot = {}
	with open(pdb_path) as pdb_f:
		for line in pdb_f:
			if not line:
				continue
			if line.startswith('DBREF'):
				fields = line.split()
				db = fields[5]
				if unp != 'UNP':
					continue
				pdbch = fields[2]
				unp_name = fields[6]
				unp_init = int(fields[8])
				unp_end = int(fields[9])
				uniprot[pdbch] = (unp_name, (unp_init, unp_end))
				continue
			if not line.startswith('ATOM'):
				continue
			res3 = line[17:20]
			ch = line[21]
			if ch not in pdb_seqs:
				pdb_seqs[ch] = ''
			res1 = from3to1(res3)
			altloc = line[16].strip()
			if not altloc:
				pdbseqs[ch] += res1
				resid = int(line[22:26].strip())
	if uniprot_spec:
		return pdb_seqs, uniprot
	return pdb_seqs
