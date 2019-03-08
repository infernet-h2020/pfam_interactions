import os
import sys
import time
import subprocess
import argparse
import pickle
import Bio.PDB
import numpy as np

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


	return conversion, uniprot_id2name, dca_resid


def download_pfam_files(pfam_acc, folder, msa_type, only_name=False, version=None):
	if version:
		print("THE version OPTION IS NOT YET IMPLEMENTED")
		exit(1)
	else:	# Assumes you need the latest version
		if msa_type == 'hmm':
			if not only_name:
				os.system("wget https://pfam.xfam.org/family/{0}/hmm -O {1}/{0}.hmm &> /dev/null".format(pfam_acc, folder))
			return folder + "/{0}.hmm".format(pfam_acc)
		else:
			if not only_name:
				os.system("wget https://pfam.xfam.org/family/{0}/alignment/{1} -O {2}/{0}_{1}.stockholm &> /dev/null".format(pfam_acc, msa_type, folder))
			return folder + "/{0}_{1}.stockholm".format(pfam_acc, msa_type)


def download_pdb(pdbname, output_folder):
	os.system("wget https://files.rcsb.org/download/{0}.pdb -O {1}/{2}.pdb &>/dev/null".format(pdbname.upper(), output_folder, pdbname.lower()))
