import os
import sys
import time
import subprocess
import argparse
import pickle
import Bio.PDB
import numpy as np

def string_is_float(s, thr_log_status='ERROR'):
	global log_filename
	this_name = string_is_float.__name__

	if not type(s) == str:
		stp_err = ('CRITICAL', this_name, "input type {0}, while str expected".format(type(s)))
		print_log(stp_err, thr_log_status=thr_log_status, log_filename=log_filename)
	try:
		float(s)
		return True
	except ValueError:
		return False


def string_is_int(s, thr_log_status='ERROR'):
	global log_filename
	this_name = string_is_int.__name__

	if not type(s) == str:
		stp_err = ('CRITICAL', this_name, "input type {0}, while str expected".format(type(s)))
		print_log(stp_err, thr_log_status=thr_log_status, log_filename=log_filename)
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


def string_decode_element(element, simple_version=False, is_immutable=False, permissive=False, thr_log_status='ERROR'):
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
			nop_err = ('CRITICAL', this_name, "could not process element {0}".format(element))
	elif permissive:
		return element
	else:
		nop_err = ('CRITICAL', this_name, "could not process immutable element {0}".format(element))

	if nop_err:
		print_log(nop_err, thr_log_status=thr_log_status, log_filename=log_filename)	


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


def backmap_pfam(target_pfam_accs, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, pfam_uniprot_stockholm_relpath, msa_type="uniprot", force_download=False):
	def print_summary(pdb_dca_resids):
		dejavu = set()
		ie = []
		kprev = ('', '')
		for k in sorted(list(pdb_dca_resids.keys())):
			c, ri = k
			prv = ''
#			print(k, pdb_dca_resids[k])
			if pdb_dca_resids[k][0][0] not in dejavu:
				if kprev[0]:
					ie[-1].append(kprev)
				ie.append([pdb_dca_resids[k][0][0], (c, ri)])
				dejavu.add(pdb_dca_resids[k][0][0])
			kprev = k
	
		if kprev[0]:
			ie[-1].append(kprev)
		for p, i, e in ie:
			print(p+":", i[0], str(i[1])+"-"+str(e[1]))


	print("Backmap: ")
	# Retrieve Pfams in the wanted PDB, and associates UniProt accession names
	pfam_in_pdb = []	# 
	text = subprocess.run(['grep', '{0}'.format(pdbname.upper()), pdb_pfam_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
	for line in text:
		if not line:
			continue
	
		fields = line.split()
		pdbc = fields[2] + '_' + fields[5]
		pfam_acc = fields[3]
		uniprot_acc = fields[4]

		if pfam_acc in [x[0] for x in target_pfam_accs]:
			pfam_in_pdb.append((pfam_acc, uniprot_acc))
	pfam_in_pdb = sorted(list(set(pfam_in_pdb)))

#	print("PFAM IN PDB", pfam_in_pdb)

	pickle_filename = "." + "".join([x+"_" for x in sorted(list(set([x[0] for x in pfam_in_pdb])))]) + "on_" + pdbname + "_" + msa_type + ".pkl"
	if os.path.exists(pickle_filename):
		bundle = pickle.load(open(pickle_filename, 'rb'))
		dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues = bundle
		print_summary(pdb_dca_resids)
		print("")
		return bundle
	
	# Conversion from UniProt to DCA resIDs
	dca_uniprot_resids = {}	# UniProt to DCA resID map. dca_uniprot_resids[(Pfam_accession, UniProt_accession)] = DCA_UniProt_map
	dca_model_length = {}	# DCA model length. dca_model_length[Pfam_accession] = DCA_length
	uniprot_restypes = {}	# UniProt resID to residue types map. uniprot_restypes[UniProt_accession] = resID_resName_map
	for pfam_acc, uniprot_acc in pfam_in_pdb:
		# Download the Pfam alignment
		pfam_uniprot_stockholm_filename = download_pfam_files(pfam_acc, pfam_uniprot_stockholm_relpath, msa_type, only_name=True)	# Only to get the correct name
		if (not os.path.exists(pfam_uniprot_stockholm_filename)) or force_download:
			pfam_uniprot_stockholm_filename = download_pfam_files(pfam_acc, pfam_uniprot_stockholm_relpath, msa_type)	# If it must be downloaded
		if not os.path.exists(pfam_uniprot_stockholm_filename):
			print("ERROR: could not find Pfam MSA")
			exit(1)
		pfam_uniprot_stockholm_filename = pfam_uniprot_stockholm_relpath + pfam_acc + '_' + msa_type + '.stockholm'
		# WARNING: This line depends on the type of Stockholm file
		text = subprocess.run(['grep', '^{0}.'.format(uniprot_acc), pfam_uniprot_stockholm_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
#		print("grep ^{0}. {1}".format(uniprot_acc, pfam_uniprot_stockholm_filename))
		if not text[0]:
			print("NOTICE: there are no sequences with the wanted uniprot access name. Request:")
			print('grep', '^{0}.'.format(uniprot_acc), pfam_uniprot_stockholm_filename)
			continue
		if not uniprot_acc in uniprot_restypes:
			uniprot_restypes[uniprot_acc] = {}
		for line in text:
			if not line:
				continue
#			print(line)
			conversion, uniprot_id2name, dca_model_length[pfam_acc] = convindex_uniprot_dca__format_stockholm(line)
			if (pfam_acc, uniprot_acc) not in dca_uniprot_resids:
				dca_uniprot_resids[(pfam_acc, uniprot_acc)] = []
			dca_uniprot_resids[(pfam_acc, uniprot_acc)].append(conversion)
			if uniprot_acc not in uniprot_restypes:
				uniprot_restypes[uniprot_acc] = uniprot_id2name
			else:
				for k in uniprot_id2name:
					uniprot_restypes[uniprot_acc][k] = uniprot_id2name[k]
	
	# Conversion from PDB to UniProt and viceversa (WARNING: they are NOT 1-to-1 correspondances)
	uniprot_pdb_resids = {}	# UniProt to PDB (each entry is a list) uniprot_pdb_resids[UniProt_accession][UniProt_resID] = [(PDB_chain, PDB_resID), ...]
	pdb_uniprot_resids = {}	# PDB to UniProt (each entry is a list) pdb_uniprot_resids[(PDB_chain, PDB_resID)] = [(UniProt_accession, UniProt_resID), ...]
	text = subprocess.run(['grep', '{0}'.format(pdbname.upper()), pdb_uniprot_res_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
	for line in text:
		if not line:
			continue
		fields = line.split()
		if fields[6] == '1':
			present = True
		else:
			present = False
		if not present:
			continue
		chain = fields[1]
#		if inch1 and chain != inch1 and chain != inch2:
#			continue
		pdb_resid = int(fields[4])
		altloc = fields[5]
		uniprot_acc = fields[8]
		uniprot_resid = int(fields[10])
		if uniprot_acc not in uniprot_pdb_resids:
			uniprot_pdb_resids[uniprot_acc] = {}
		if uniprot_resid not in uniprot_pdb_resids[uniprot_acc]:
			uniprot_pdb_resids[uniprot_acc][uniprot_resid] = set()
		uniprot_pdb_resids[uniprot_acc][uniprot_resid].add((chain, pdb_resid))
		if (chain, pdb_resid) not in pdb_uniprot_resids:
			pdb_uniprot_resids[(chain, pdb_resid)] = []
		pdb_uniprot_resids[(chain, pdb_resid)].append((uniprot_acc, uniprot_resid))
	
	# Conversion DCA to PDB and viceversa
	dca_pdb_resids = {}	# DCA to PDB. dca_pdb_resids[(Pfam_accession, UniProt_accession)] = [[[DCA_resid, PDB_resid], ...], ]
	pdb_dca_resids = {}	# PDB to DCA. pdb_dca_resids[()]
	chain_multiplicity = {}	# Times a Pfam is contained in each PDB chain. chain_multiplicity[Pfam_accession][PDB_chain] = multiplicity
	for pfam_acc, uniprot_acc in pfam_in_pdb:
		new_mainlist = []
		if pfam_acc not in chain_multiplicity:
			chain_multiplicity[pfam_acc] = {}
		for convlist in dca_uniprot_resids[(pfam_acc, uniprot_acc)]:	# WARNING: The first entry contains the initial and final UniProt resIDs
			new_convlist = []
			chains_present = set()
			for i, p in enumerate(convlist):
				if i == 0:
					initial_uniprot_resid, final_uniprot_resid = p
					continue
				dca_resid, uniprot_resid = p
				if uniprot_resid in uniprot_pdb_resids[uniprot_acc]:
					pdb_resid = uniprot_pdb_resids[uniprot_acc][uniprot_resid]
					new_convlist.append([dca_resid, pdb_resid])
#					print(dca_resid, uniprot_resid, pdb_resid)
					for c, ri in pdb_resid:
						if c not in chains_present:
							if c not in chain_multiplicity[pfam_acc]:
								chain_multiplicity[pfam_acc][c] = 0
							chain_multiplicity[pfam_acc][c] += 1
							chains_present.add(c)
#							print(pfam_acc, c, chain_multiplicity[pfam_acc][c])
						if (c, ri) not in pdb_dca_resids:
							pdb_dca_resids[(c, ri)] = []
						
						el = (pfam_acc + '_' + c + str(chain_multiplicity[pfam_acc][c]), dca_resid)
						if el not in pdb_dca_resids[(c, ri)]:
							pdb_dca_resids[(c, ri)].append((pfam_acc + '_' + c + str(chain_multiplicity[pfam_acc][c]), dca_resid))
			new_mainlist.append(new_convlist)
		dca_pdb_resids[(pfam_acc, uniprot_acc)] = new_mainlist[:]
	
	# List the residues with coordinates in each PDB chain
	allowed_residues = {}
	for x in sorted(list(pdb_dca_resids.keys())):
		chain = x[0]
		resid = x[1]
		if chain not in allowed_residues:
			allowed_residues[chain] = set()
		allowed_residues[chain].add(resid)

#	print(pickle_filename)
	pickle.dump((dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues), open(pickle_filename, 'wb'))
	print_summary(pdb_dca_resids)	

	print("")
	return dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues 


def compute_accessibility(pdbname, pdb_path, pdb_dca_resids, allowed_residues, vmd_path, results_folder):
	print("SASA:")
	residues_linear = []
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbname, pdb_path)
	distance_threshold = 20
	for model in structure:
		for chain1 in allowed_residues:
			c1 = model[chain1]
			for r1 in c1:
				resid1 = r1.id[1]
				if resid1 not in allowed_residues[chain1]:
					continue
				residues_linear.append((chain1, resid1))

	N = len(residues_linear)

#	print(residues_linear)
#	print(allowed_residues)

	SASA = {}
	for i1 in range(N):
		c1, r1 = residues_linear[i1]
		for pfam_acc_ann1, dca_i1 in pdb_dca_resids[(c1, r1)]:
#			if dca_i1 in accessibility_dca_resids:
			pfam_acc1, cmult1 = pfam_acc_ann1.split('_')
			os.system("sed 's/XXX/chain {0} and resid {1}/' fixed_sasa_template.tcl > fixed_sasa.tcl".format(c1, r1))
#			subprocess.run(["/Applications/VMD\ 1.9.3.app/Contents/vmd/vmd_MACOSXX86", pdb_files_ext_path + pdbname.lower() + '.pdb', "-dispdev", "text", "-e", "fixed_sasa.tcl"], shell=True)
			os.system(vmd_path + " " + pdb_path + " -dispdev text -e fixed_sasa.tcl &> /dev/null")
			with open("SASA_output.dat") as output_file:
				for line in output_file:
					if not dca_i1 in SASA:
						SASA[dca_i1] = []
					SASA[dca_i1].append((pfam_acc_ann1, pdbname+' '+c1+' '+str(r1), float(line.strip())))
					break
	os.remove("SASA_output.dat")

	acessibility_filename = results_folder + "{0}_SASA.txt".format(pdbname)
	with open(acessibility_filename, 'w') as accessibility_file:
		for dca_i in SASA:
			for x in SASA[dca_i]:
				accessibility_file.write("{0}\t{1}\t{2}\t{3}\n".format(dca_i, x[0], x[1], x[2]))
	print(acessibility_filename)
	return SASA
	

def compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder):
	print("Interactions:")
	t = time.time()
	bb_names = {'N', 'O', 'OXT', 'C'}
	distance_threshold = 20
	print("\tcalculating distances\t", end='', flush=True)
	pickle_filename = "." + pdbname + "_distances.pkl"
	if os.path.exists(pickle_filename):
		N, residues_linear, distance, sc_distance, CA_distance = pickle.load(open(pickle_filename, 'rb'))
		print("(pickled!)\t", end='', flush=True)
	else:
		residues_linear = []
		distance = []
		sc_distance = []
		CA_distance = []
		parser = Bio.PDB.PDBParser(QUIET=True)
		structure = parser.get_structure(pdbname, pdb_path)
		
		for model in structure:
			for chain1 in allowed_residues:
				c1 = model[chain1]
				for r1 in c1:
					resid1 = r1.id[1]
					if resid1 not in allowed_residues[chain1]:
						continue
					residues_linear.append((chain1, resid1))
					for chain2 in allowed_residues:
						c2 = model[chain2]
						for r2 in c2:
							resid2 = r2.id[1]
							if resid2 not in allowed_residues[chain2]:
								continue
							if (chain1 == chain2 and resid1 == resid2) or (inch1 and ({chain1, chain2} != {inch1, inch2})):
								distance.append(0)
								sc_distance.append(0)
								CA_distance.append(0)
								continue
							CA_CA_d = np.linalg.norm(r1['CA'].get_vector() - r2['CA'].get_vector())
							CA_distance.append(CA_CA_d)
							if CA_CA_d > distance_threshold:
								sc_distance.append(CA_CA_d)
								distance.append(CA_CA_d)
								continue
							mindist = 100000
							mindist_sc = 100000
							for a1 in r1:
								for a2 in r2:
									d = np.linalg.norm(a1.get_vector() - a2.get_vector())
									if d < mindist:
										mindist = d
									if d < mindist_sc and (a1.id not in bb_names) and (a2.id not in bb_names):
										mindist_sc = d
							distance.append(mindist)
							sc_distance.append(mindist_sc)
			break	# consider only first model
		
		N = len(residues_linear)
		
		distance = np.array(distance).reshape((N, N))
		sc_distance = np.array(sc_distance).reshape((N, N))
		CA_distance = np.array(CA_distance).reshape((N, N))

		pickle.dump((N, residues_linear, distance, sc_distance, CA_distance), open(pickle_filename, 'wb'))

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))
	t = tf
	
	interactions, sc_interactions, CA_interactions= {}, {}, {}
	ausilia = {}
	print("\tpreparing support data for calculating interactions\t", end='', flush=True)
	for i1 in range(N):
		c1, r1 = residues_linear[i1]
		for i2 in range(N):
			c2, r2 = residues_linear[i2]
			if inch1 and {c1, c2} != {inch1, inch2}:
				continue
			for pfam_acc_ann1, dca_i1 in pdb_dca_resids[(c1, r1)]:
				pfam_acc1, cmult1 = pfam_acc_ann1.split('_')
				L1 = dca_model_length[pfam_acc1]
				for pfam_acc_ann2, dca_i2 in pdb_dca_resids[(c2, r2)]:
					pfam_acc2, cmult2 = pfam_acc_ann2.split('_')
					L2 = dca_model_length[pfam_acc2]
					if not (pfam_acc1 == pfam_in_pdb[0][0] and pfam_acc2 == pfam_in_pdb[1][0]):	# Maintains the order of Pfams chosen at the beginning (in the command line!)
						continue
					if (pfam_acc1, pfam_acc2) not in interactions:
						interactions[(pfam_acc1, pfam_acc2)], sc_interactions[(pfam_acc1, pfam_acc2)], CA_interactions[(pfam_acc1, pfam_acc2)] = {}, {}, {}
						ausilia[(pfam_acc1, pfam_acc2)] = {}
					if (cmult1, cmult2) not in interactions[(pfam_acc1, pfam_acc2)]:
						interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)], sc_interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)], CA_interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)] = -np.ones((L1,L2)), -np.ones((L1,L2)), -np.ones((L1,L2))
						ausilia[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)] = {}
					if dca_i1-1 not in ausilia[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)]:
						ausilia[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][dca_i1-1] = {}
					d, dSC, dCA = distance[i1][i2], sc_distance[i1][i2], CA_distance[i1][i2]
					### Filter
	#				if c1 == c2 and abs(r1 - r2) < 4:
	#					d = -1
					###
					interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][dca_i1-1][dca_i2-1] = d
					sc_interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][dca_i1-1][dca_i2-1] = dSC
					CA_interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][dca_i1-1][dca_i2-1] = dCA
					ausilia[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][dca_i1-1][dca_i2-1] = (c1, r1, c2, r2)

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))
	t = tf

	multiplicities = {}
	print("\tcalculating interactions\t", end='', flush=True)
	interaction_filenames = set()
	interaction_summary = []
	for x in interactions:
		for imm, mm in enumerate(list(interactions[x].keys())):
			cmult1, cmult2 = mm
			c1, c2 = cmult1[0], cmult2[0]
			if inch1 and {c1, c2} != {inch1, inch2}:
				continue
			pfam_acc1, pfam_acc2 = x
			interactions_filename = results_folder + "{0}_{1}_{2}_{3}_{4}.txt".format(pfam_acc1, pfam_acc2, pdbname, cmult1, cmult2)
			if pfam_acc1 == pfam_acc2 and cmult1 == cmult2:
				attr = 'same'
				offset = 0
			elif c1 == c2:
				attr = 'intra'
				offset = dca_model_length[pfam_acc1]
			else:
				attr = 'inter'
				offset = dca_model_length[pfam_acc1]
			with open(interactions_filename, 'w') as interactions_file:
				L1, L2 = interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)].shape
				for i in range(L1):
					for j in range(L2):
						val, sc_val, CA_val = interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][i][j], sc_interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][i][j], CA_interactions[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][i][j]
						if val == -1:
							interactions_file.write("{0:6}\t{1:6}\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\n".format(i+1, j+1+offset))
						else:
							ac1, ar1, ac2, ar2 = ausilia[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][i][j]
							try:
								uniprot_resids1 = "".join([x[0] + "_" + str(x[1]) + "_" + uniprot_restypes[x[0]][x[1]] + "," for x in pdb_uniprot_resids[(ac1, ar1)]])[:-1]
							except:
								uniprot_resids1 = "None"
							try:
								uniprot_resids2 = "".join([x[0] + "_" + str(x[1]) + "_" + uniprot_restypes[x[0]][x[1]] + "," for x in pdb_uniprot_resids[(ac2, ar2)]])[:-1]
							except:
								uniprot_resids2 = "None"
	
							interactions_file.write("{0:6}\t{1:6}\t{2:4}\t{3:4}\t{4:4}\t{5:4}\t{6:10.4}\t{7:10.4}\t{8:10.4}\t{9}\t{10}\n".format(i+1, j+1+offset, ac1, ar1, ac2, ar2, val, sc_val, CA_val, uniprot_resids1, uniprot_resids2))
			if (pfam_acc1, pfam_acc2) not in multiplicities:
				multiplicities[(pfam_acc1, pfam_acc2)] = []
			if (cmult1, cmult2) not in multiplicities[(pfam_acc1, pfam_acc2)]:
				multiplicities[(pfam_acc1, pfam_acc2)].append((cmult1, cmult2))
			multinum = multiplicities[(pfam_acc1, pfam_acc2)].index((cmult1, cmult2))
			interaction_summary.append((pfam_acc1, pfam_acc2, multinum, interactions_filename, attr))
			interaction_filenames.add(interactions_filename)

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))

	for pfam_acc1, pfam_acc2, multinum, interactions_filename, attr in interaction_summary:
		print(pfam_acc1, pfam_acc2, multinum, interactions_filename, attr)

	if not interaction_filenames:
		print("")
	print("")
	return interaction_filenames


def calculate_matches(options, pdbname=''):
	if not pdbname:
		pdbname = options['pdbname']
	pdb_pfam_filename = options['pdb_pfam_filename']
	inpfam = options['inpfam']
	inpfam1 = options['inpfam1']
	inpfam2 = options['inpfam2']

	print("Matches:")
	pfam_in_pdb = ['', '']
	text = subprocess.run(['grep', '{0}'.format(pdbname.upper()), pdb_pfam_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
#	print("grep {0} {1}".format(pdbname, pdb_pfam_filename))
#	print(text)
	for line in text:
		if not line:
			continue
		fields = line.split()
		pdbc = fields[2] + '_' + fields[5]
		pfam_acc = fields[3]
		uniprot_acc = fields[4]
		if (inpfam and pfam_acc != inpfam) or (inpfam1 and (pfam_acc != inpfam1 and pfam_acc != inpfam2)):
			continue
		elif inpfam1 and pfam_acc == inpfam1:
			pfam_in_pdb[0] = (pfam_acc, uniprot_acc)
		elif inpfam1 and pfam_acc == inpfam2:
			pfam_in_pdb[1] = (pfam_acc, uniprot_acc)
		elif inpfam:
			pfam_in_pdb.append((pfam_acc, uniprot_acc))
		print(line)
#		print(pfam_in_pdb)
	pfam_in_pdb_correct = []
	for p in pfam_in_pdb:
		if p:
			pfam_in_pdb_correct.append(p)
	pfam_in_pdb = pfam_in_pdb_correct
	print("")

	# Consistency checks	
	if (inpfam1 and (inpfam1 not in [x[0] for x in pfam_in_pdb] or inpfam2 not in [x[0] for x in pfam_in_pdb])) or inpfam and (inpfam not in [x[0] for x in pfam_in_pdb]):
		print("ERROR: One or more query Pfams was not found in target PDB")
		exit(1)
	if not pfam_in_pdb:
		print("No Pfams selected!")
		exit(1)

	return pfam_in_pdb


def main_parser():
	options = {}	

	# Hardcoded paths
	options['database_files_relpath'] = "../examples/Pfam_32_uniprot/database_files/"	# Folder where the external files from the Pfam database are kept
	options['pdb_pfam_filename'] = options['database_files_relpath'] + "pdb_pfamA_reg.txt"	# Register cennecting PDB and Pfam information
	options['pdb_uniprot_res_filename'] = options['database_files_relpath'] + "pdb_residue_data.txt"	# Register with residue-by-residue details
	options['pfam_uniprot_filename'] = options['database_files_relpath'] + "uniprot_reg_full.txt"	# Register connecting UniProt and Pfam information	
	options['pfam_uniprot_stockholm_relpath'] = "../examples/Pfam_32_uniprot/alns_stockholm/"	# Folder where alignments are stored
	options['results_folder'] = "../examples/results/"	# Folder where results are put
#	pdb_files_ext_path = "/media/sarti/data/Work/databases/PDB/"
	options['pdb_files_ext_path'] = "../examples/PDB/"
	options['pfam_pfam_filename'] = "../complete_database/PfamPfam_contacts_through_all_PDBs_table.txt"

	# Option parser
	parser = argparse.ArgumentParser()

	parser.add_argument('-pdb', '--pdbname', nargs='?')
	parser.add_argument('-a', '--accessibility', action='store_const', const='True', default='False')
	parser.add_argument('-ip', '--inpfam', nargs='?')
	parser.add_argument('-ip1', '--inpfam1', nargs='?')
	parser.add_argument('-ip2', '--inpfam2', nargs='?')
	parser.add_argument('-ic1', '--inch1', nargs='?')
	parser.add_argument('-ic2', '--inch2', nargs='?')
	parser.add_argument('-ptype', '--pfam_msa_type', nargs='?')
	parser.add_argument('-vmd', '--vmd_path', nargs='?')
	parser.add_argument('-fd', '--force_download', nargs='?')
	parser.add_argument('-mind', '--min_dist', nargs='?')
	parser.add_argument('-dca', '--dca_filename', nargs='?')
	parser.add_argument('-intra', '--only_intra', action='store_const', const='True', default='False')
	parser.add_argument('-inter', '--only_inter', action='store_const', const='True', default='False')

	# Default values for optional arguments
	parser.set_defaults(pdbname = 'None')
	parser.set_defaults(accessibility = 'False')
	parser.set_defaults(inpfam = 'None')
	parser.set_defaults(inpfam1 = 'None')
	parser.set_defaults(inpfam2 = 'None')
	parser.set_defaults(inch1 = 'None')
	parser.set_defaults(inch2 = 'None')
	parser.set_defaults(pfam_msa_type = 'uniprot')
	parser.set_defaults(vmd_path = '/Applications/VMD\ 1.9.3.app/Contents/vmd/vmd_MACOSXX86')
	parser.set_defaults(force_download = 'False')
	parser.set_defaults(min_dist = 'None')
	parser.set_defaults(dca_filename = 'None')

	# Parse options
	parsed = parser.parse_args()
#	print(parsed)
	for x in parsed.__dict__:
#		if x in parsed.__dict__ and not parsed.__dict__[x]:
#			print(x)
#			continue
#		print(x, parsed.__dict__[x], type(parsed.__dict__[x]))
		options[x] = string_decode_element(parsed.__dict__[x], permissive=True, simple_version=True)
#	print(options)
#	print(parsed)
#	print(parsed.__dict__)

	pfam_msa_types = {'uniprot' : 'uniprot', 'full' : 'full', 'rp75' : 'rp75'}
	if options['pfam_msa_type'] not in pfam_msa_types:
		print("ERROR: Pfam MSA type " + options['pfam_msa_type'] + " not defined")
		exit(1)
	options['msa_type'] = pfam_msa_types[options['pfam_msa_type']]

	if type(options['inpfam']) == type(None) and type(options['inpfam1']) == type(None) and type(options['inpfam2']) == type(None):
		print("ERROR: Specify at least one Pfam accession code")
		exit(1)

	for k in ['inpfam', 'inpfam1', 'inpfam2']:
		x = options[k]
		if type(x) != type(None) and (len(x) != 7 or x[:2] != "PF" or (not string_is_int(x[2:]))):
			print("ERROR: argument of option {0} is badly formatted".format(k))
			exit(1)

	if options['accessibility'] == True:
		main_accessibility(options)
	elif type(options['min_dist']) != type(None):
		main_mindist(options)
	else:
		main_interactions(options)


def main_interactions(options):
	pdbname = options['pdbname']
	pdb_pfam_filename = options['pdb_pfam_filename']
	pdb_uniprot_res_filename = options['pdb_uniprot_res_filename']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	inch1 = options['inch1']
	inch2 = options['inch2']
	results_folder = options['results_folder']

	# Download PDB if needed
	pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
	if not os.path.exists(pdb_path):
		if not os.path.exists(options['pdb_files_ext_path']):
			print("ERROR: PDB path not found")
			exit(1)
		download_pdb(pdbname, options['pdb_files_ext_path'])
	if not os.path.exists(pdb_path):
		print("ERROR: PDB " + pdbname + " not found")
		exit(1)
	options['pdb_path'] = pdb_path

	pfam_in_pdb = calculate_matches(options)

	dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues = backmap_pfam(pfam_in_pdb, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, pfam_uniprot_stockholm_relpath, msa_type=msa_type, force_download=force_download)

	compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder)


def download_pdb(pdbname, output_folder):
	os.system("wget https://files.rcsb.org/download/{0}.pdb -O {1}/{2}.pdb &>/dev/null".format(pdbname.upper(), output_folder, pdbname.lower()))
	

def main_accessibility(options):
	pdbname = options['pdbname']
	pdb_pfam_filename = options['pdb_pfam_filename']
	pdb_uniprot_res_filename = options['pdb_uniprot_res_filename']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	vmd_path = options['vmd_path']
	results_folder = options['results_folder']
	
	if (not vmd_path):
		print("ERROR: SASA calculation needs VMD to run. Specify a path for VMD")
		exit(1)

	# Download PDB if needed
	pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
	if not os.path.exists(pdb_path):
		if not os.path.exists(options['pdb_files_ext_path']):
			print("ERROR: PDB path not found")
			exit(1)
		download_pdb(pdbname, options['pdb_files_ext_path'])
	if not os.path.exists(pdb_path):
		print("ERROR: PDB " + pdbname + " not found")
		exit(1)
	options['pdb_path'] = pdb_path

	pfam_in_pdb = calculate_matches(options)
	dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues = backmap_pfam(pfam_in_pdb, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, pfam_uniprot_stockholm_relpath, msa_type=msa_type, force_download=force_download)

	SASA = compute_accessibility(pdbname, pdb_path, pdb_dca_resids, allowed_residues, vmd_path, results_folder)


def main_mindist(options):
	mindist = options['min_dist']
	inpfam1 = options['inpfam1']
	inpfam2 = options['inpfam2']
	pfam_pfam_filename = options['pfam_pfam_filename']
	pdb_pfam_filename = options['pdb_pfam_filename']
	pdb_uniprot_res_filename = options['pdb_uniprot_res_filename']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	pdb_files_ext_path = options['pdb_files_ext_path']
	results_folder = options['results_folder']
	dca_filename = options['dca_filename']
	only_intra = options['only_intra']
	only_inter = options['only_inter']
	inch1 = ''
	inch2 = ''

	if not dca_filename:
		print("ERROR: mindist needs a precomputed plmdca filename")
	
	mind_pdbs = []
	if mindist == 'all':
		text = subprocess.run(["grep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam1, pfam_pfam_filename, inpfam2)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
		for line in text:
			if not line:
				continue
			if line.strip() not in mind_pdbs:
				mind_pdbs.append(line.strip())
		mind_pdbs = sorted(list(mind_pdbs))
	else:
		if os.path.exists(mindist):
			with open(mindist) as mindist_file:
				for line in mindist_file:
					if not line:
						continue
					if line.strip() not in mind_pdbs:
						mind_pdbs.append(line.strip())
			mind_pdbs = sorted(list(mind_pdbs))
		else:
			print("ERROR: The file specifying the PDB names for the mindist analysis cannot be found")
			exit(1)

	print("\nConsidering the following PDBs:")
	print("".join([x+", " for x in mind_pdbs])[:-2]+"\n")

	interaction_filenames = set()
	failed_pdbs = set()
	for pdbname in mind_pdbs:
		print("PDB: ", pdbname, '-'*150)
		pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
		try:
			pfam_in_pdb = calculate_matches(options, pdbname=pdbname)
		except:
			failed_pdbs.add(pdbname)
			continue

		# Download PDB if needed
		pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
		if not os.path.exists(pdb_path):
			if not os.path.exists(options['pdb_files_ext_path']):
				print("ERROR: PDB path not found")
				exit(1)
			download_pdb(pdbname, options['pdb_files_ext_path'])
		if not os.path.exists(pdb_path):
			print("ERROR: PDB " + pdbname + " not found")
			exit(1)

		dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues = backmap_pfam(pfam_in_pdb, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, pfam_uniprot_stockholm_relpath, msa_type=msa_type, force_download=force_download)

		int_filenames = compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder)	
		interaction_filenames |= int_filenames	

	new_mind_pdbs = []
	for x in mind_pdbs:
		if x in failed_pdbs:
			continue
		new_mind_pdbs.append(x)
	mind_pdbs = new_mind_pdbs

	text = []
	for pdbname in mind_pdbs:
		text += subprocess.run(["grep {0} {1} | grep {2} | grep {3} | awk '{{if ($1==\"{0}\") {{print \"{3}\", $4, $5}} else {{print \"{3}\", $5, $4}} }}'".format(inpfam1, pfam_pfam_filename, inpfam2, pdbname)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
#		subprocess.run(["grep {0} {1} | grep {2} | grep {3}".format(inpfam1, pfam_pfam_filename, inpfam2, pdbname)], shell=True)

	pairs = []
	for line in text:
		if not line:
			continue
		pdb, p1, p2 = line.strip().split()
		if only_intra and p1 == p2:
			pairs.append((pdb, p1, p2))
		elif only_inter and p1 != p2:
			pairs.append((pdb, p1, p2))
		elif not (only_intra or only_inter):
			pairs.append((pdb, p1, p2))

	filenames_list = []
	for pair in pairs:
		fname = results_folder + "/{0}_{1}_{2}_{3}*_{4}*.txt".format(inpfam1, inpfam2, pair[0].lower(), pair[1], pair[2])
		text = subprocess.run(["ls {0} 2>/dev/null".format(fname)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
		for line in text:
			if not line:
				continue
			filenames_list.append(line.strip())
	
#	print("FILES", filenames_list)

	with open(filenames_list[0]) as one_file:
		for line in one_file:
			if not line:
				continue
			offset = int(line.split()[0])

	imax = 50
	i = 0
	out_filename = results_folder + "/{0}_{1}_correlation.txt".format(inpfam1, inpfam2)
	if only_intra:
		warner = " (only intra-chain)"
	elif only_inter:
		warner = " (only inter-chain)"
	else:
		warner = ""
	print("\n\nStrongest DCA signals mapped on closest residue pairs{0}:".format(warner))
	with open(dca_filename, 'r') as dca_file:
		with open(out_filename, 'w') as out_file:
			for nl, line in enumerate(dca_file):
				if not line:
					continue
				fields = line.split()
				if not (int(fields[0]) <= offset and int(fields[1]) > offset):
					continue
				outs = []
				for res_filename in filenames_list:
					text = subprocess.run(["grep '^\s*{0}\s\s*{1}\s' {2}".format(fields[0], fields[1], res_filename)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
					if text:
						text = text[:-1]
						outs.append((res_filename, text, float(fields[2])))
				mind = 1000000
				argmind = ""
				for o in outs:
					d = o[1].split()[6]
					if d == "None":
						continue
					elif float(d) < mind:
						mind = float(d)
						argmind = o
				if argmind:
					i += 1
					print("{0}\t{1}\t{2:.5f}\t{3}".format(argmind[0], argmind[1], argmind[2], nl))
					out_file.write("{0}\t{1}\t{2:.5f}\t{3}\n".format(argmind[0], argmind[1], argmind[2], nl))
					if i >= imax:
						break


if __name__ == "__main__":
	main_parser()
