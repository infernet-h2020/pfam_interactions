from support import *

def backmap_pfam(target_pfam_accs, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, pfam_uniprot_stockholm_relpath, cache_folder, msa_type="uniprot", force_download=False):
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

	pickle_filename = cache_folder + "." + "".join([x+"_" for x in sorted(list(set([x[0] for x in pfam_in_pdb])))]) + "on_" + pdbname + "_" + msa_type + ".pkl"
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