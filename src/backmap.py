from support import *

def backmap_pfam(target_pfam_accs, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type="uniprot", force_download=False, complete_backmap=False):
	def print_summary(pdb_dca_resids, pdb_uniprot_resids, backmap_filename):
		dejavu = set()
		ie = []
		kprev = ('', '')

		with open(backmap_filename, 'w') as bkmp_file:
			for k in sorted(list(pdb_dca_resids.keys())):
				c, ri = k
				prv = ''
	
				for inst in pdb_dca_resids[k]:
					bkmp_file.write("{0}\t{1:5d}\t{2}\t{3:5d}\n".format(c, ri, inst[0], inst[1]))
	
				if pdb_dca_resids[k][0][0] not in dejavu:
					if kprev[0]:
						ie[-1].append(kprev)
					ie.append([pdb_dca_resids[k][0][0], (c, ri)])
					dejavu.add(pdb_dca_resids[k][0][0])
				kprev = k
	
		if kprev[0]:
			ie[-1].append(kprev)

		backmap_table = []
		for p, i, e in ie:
			up_str = ""
			up = []
			for ipur in range(len(pdb_uniprot_resids[(i[0], i[1])])):
				if up_str:
					up_str += " "
				up_str += pdb_uniprot_resids[(i[0], i[1])][ipur][0] + "_" + str(pdb_uniprot_resids[(i[0], i[1])][ipur][1]) + "-" + str(pdb_uniprot_resids[(e[0], e[1])][ipur][1])
				up.append((pdb_uniprot_resids[(i[0], i[1])][ipur][0], (pdb_uniprot_resids[(i[0], i[1])][ipur][1], pdb_uniprot_resids[(e[0], e[1])][ipur][1])))
			print(p+":", i[0], str(i[1])+"-"+str(e[1]), up_str)
			backmap_table.append([p, (i[0], (i[1], e[1])), tuple(up)])	# pfam_acc_ann, chain, init res, end res, mapping uniprot

		print("Find the complete backmapping in {0}".format(backmap_filename))
#		print(backmap_table)
		return backmap_table

	def vmd_backmap_script(backmap_table, tcl_output_filename, pdb_path):
		colors = [0,1,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]	# no grey (2), no white (8)
		with open(tcl_output_filename, 'w') as tcl_output_file:
			tcl_output_file.write("display resetview\n")
			tcl_output_file.write("display projection Orthographic\n")	# Ortographic projection
			tcl_output_file.write("display depthcue off\n")			# No depth shade
			tcl_output_file.write("color Display Background white\n")
			tcl_output_file.write("mol new {{{0}}} type {{pdb}} first 0 last -1 step 1 waitfor 1\n".format(pdb_path))
			lastipfd = 0
			for ipfd, pfd in enumerate(backmap_table):
				tcl_output_file.write("mol addrep 0\n")
				chname = pfd[1][0]
				init = int(pfd[1][1][0])
				end = int(pfd[1][1][1])
				tcl_output_file.write("# Pfam {0}\n".format(pfd[0]))
				tcl_output_file.write("mol modselect {0} 0 chain {1} and resid {2} to {3}\n".format(ipfd, chname, init, end))
				tcl_output_file.write("mol modstyle {0} 0 Licorice 0.300000 12.000000 12.000000\n".format(ipfd))
				tcl_output_file.write("mol modcolor {0} 0 ColorID {1}\n".format(ipfd, colors[ipfd]))
				lastipfd = ipfd
			tcl_output_file.write("# Trace \n")
			tcl_output_file.write("mol addrep 0\n")
			tcl_output_file.write("mol modstyle {0} 0 Trace 0.300000 12.000000\n".format(lastipfd+1))
			tcl_output_file.write("mol modcolor {0} 0 ColorID 2\n".format(lastipfd+1))

	# Read index
	pdb_uniprot_index = {}
	src_path = os.path.dirname(os.path.abspath(__file__)) + '/'
	with open(pdb_uniprot_res_index_filename) as index_file:
		for line in index_file:
			if not line.strip():
				continue
			fields = line.split()
			#pdb_uniprot_index[fields[0]] = (src_path + fields[1], fields[2])
			pdb_uniprot_index[fields[0]] = (fields[1], fields[2])


	print("Backmap: ")
	# Retrieve Pfams in the wanted PDB, and associates UniProt accession names
	pfam_in_pdb = []	# 
#	print('grep {0} {1}'.format(pdbname.upper(), pdb_pfam_filename))
	text = subprocess.run(['zgrep', '{0}'.format(pdbname.upper()), pdb_pfam_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
	for line in text:
		if not line:
			continue
		fields = line.split()
		pdb = fields[2]
		if pdb != pdbname.upper():
			continue
		pdbc = fields[2] + '_' + fields[5]
		pfam_acc = fields[3]
		uniprot_acc = fields[4]
		if complete_backmap:
			pfam_in_pdb.append((pfam_acc, uniprot_acc))
		elif pfam_acc in [x[0] for x in target_pfam_accs]:
			pfam_in_pdb.append((pfam_acc, uniprot_acc))
	pfam_in_pdb = sorted(list(set(pfam_in_pdb)))

#	print("PFAM IN PDB", pfam_in_pdb)

	if complete_backmap:
		pickle_filename = cache_folder + ".complete_backmap_of_" + pdbname + "_" + msa_type + '_v' + str(version) + ".pkl"
		backmap_filename = results_folder + "complete_backmap_of_" + pdbname + "_" + msa_type + '_v' + str(version) + ".txt"
		
	else:
		pickle_filename = cache_folder + "." + "".join([x+"_" for x in sorted(list(set([x[0] for x in pfam_in_pdb])))]) + "on_" + pdbname + "_" + msa_type + '_v' + str(version) + ".pkl"
		backmap_filename = results_folder + "".join([x+"_" for x in sorted(list(set([x[0] for x in pfam_in_pdb])))]) + "on_" + pdbname + "_" + msa_type + '_v' + str(version) + ".txt"
	
#	print(pickle_filename, backmap_filename)
	if os.path.exists(pickle_filename):
		bundle = pickle.load(open(pickle_filename, 'rb'))
		dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle
		backmap_table = print_summary(pdb_dca_resids, pdb_uniprot_resids, backmap_filename)
		if len(target_pfam_accs) == 1:
			tcl_output_filename = results_folder + "/{0}_domain_visualization.tcl".format(target_pfam_accs[0][0])
		elif len(target_pfam_accs) == 2:
			tcl_output_filename = results_folder + "/{0}_{1}_domain_visualization.tcl".format(target_pfam_accs[0][0], target_pfam_accs[1][0])
		vmd_backmap_script(backmap_table, tcl_output_filename, pdb_path)
		print("To visualize the domains: vmd -e {0}".format(tcl_output_filename))
		print("")
		return bundle

	print("\tUniProt resIDs >> DCA resIDs\t", end='', flush=True)
	t = time.time()
	
	# Conversion from UniProt to DCA resIDs
	dca_uniprot_resids = {}	# UniProt to DCA resID map. dca_uniprot_resids[(Pfam_accession, UniProt_accession)] = DCA_UniProt_map
	dca_model_length = {}	# DCA model length. dca_model_length[Pfam_accession] = DCA_length
	uniprot_restypes = {}	# UniProt resID to residue types map. uniprot_restypes[UniProt_accession] = resID_resName_map
	delete_uniprot_acc = set()
	for pfam_acc, uniprot_acc in pfam_in_pdb:
		# Download the Pfam alignment
		pfam_uniprot_stockholm_filename = download_pfam_files(pfam_acc, pfam_uniprot_stockholm_relpath, msa_type, version, only_name=True)	# Only to get the correct name

#		if (not os.path.exists(pfam_uniprot_stockholm_filename)) or force_download:
#			pfam_uniprot_stockholm_filename = download_pfam_files(pfam_acc, pfam_uniprot_stockholm_relpath, msa_type, version)	# If it must be downloaded
		if not os.path.exists(pfam_uniprot_stockholm_filename):
			print(pfam_uniprot_stockholm_filename)
			print("\nERROR: could not find Pfam MSA")
			exit(1)
		# WARNING: This line depends on the type of Stockholm file
		text = subprocess.run(['zgrep', '^{0}.'.format(uniprot_acc), pfam_uniprot_stockholm_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
#		print("grep ^{0}. {1}".format(uniprot_acc, pfam_uniprot_stockholm_filename))
		if not text[0]:
			print("\nNOTICE: there are no sequences with the wanted uniprot access name. Request:")
			print('zgrep', '^{0}.'.format(uniprot_acc), pfam_uniprot_stockholm_filename)
			delete_uniprot_acc.add(uniprot_acc)
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

	new_pfam_in_pdb = []
	for pfam_acc, uniprot_acc in pfam_in_pdb:
		if uniprot_acc not in delete_uniprot_acc:
			new_pfam_in_pdb.append((pfam_acc, uniprot_acc))
	pfam_in_pdb = new_pfam_in_pdb[:]

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))
	print("\tPDB <<>> UniProt\t", end='', flush=True)
	t = tf

	# Conversion from PDB to UniProt and viceversa (WARNING: they are NOT 1-to-1 correspondances)
	uniprot_pdb_resids = {}	# UniProt to PDB (each entry is a list) uniprot_pdb_resids[UniProt_accession][UniProt_resID] = [(PDB_chain, PDB_resID), ...]
	pdb_uniprot_resids = {}	# PDB to UniProt (each entry is a list) pdb_uniprot_resids[(PDB_chain, PDB_resID)] = [(UniProt_accession, UniProt_resID), ...]
#	print('grep', '{0}'.format(pdbname.upper()), pdb_uniprot_res_filename)
#	text = subprocess.run(['grep', '{0}'.format(pdbname.upper()), pdb_uniprot_res_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')

	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbname, pdb_path)
	if not structure:
		print("BIOPYTHON PARSER DID NOT PARSE")
		return None

#	print("PDB", pdb_path, pdbname)

	valid_residues = {}
#	print("VALID RESIDUES", pdbname)
	for chain in structure[0]:
		valid_residues[chain.id] = [x.id[1] for x in chain if (not x.id[0].strip()) and not( x.id[2].strip()) and 'CA' in x]
#		print(pdbname, chain.id, valid_residues[chain.id])

	indexed_pdb_uniprot_res_filename, ind = pdb_uniprot_index[pdbname.upper()]
	with open(indexed_pdb_uniprot_res_filename) as indexed_pdb_uniprot_res_file:
		for ln, line in enumerate(indexed_pdb_uniprot_res_file):
			if not line or ln < int(ind)-1:
				continue
			fields = line.split()
			if fields[0] != pdbname.upper():
				continue
			if fields[6] == '1':
				present = True
			else:
				present = False
			if not present:
#				print("Not present")
				continue
			chain = fields[1]
			if chain not in valid_residues:
#				print("chain not valid")
				continue
#			if inch1 and chain != inch1 and chain != inch2:
#				continue
			pdb_resid = int(fields[4])
			if pdb_resid not in valid_residues[chain]:
#				print("resid not valid")
				continue
			altloc = fields[5]
			uniprot_acc = fields[8]
			uniprot_resid = int(fields[10])
#			print("BOEUF", uniprot_acc)
			if uniprot_acc not in uniprot_pdb_resids:
				uniprot_pdb_resids[uniprot_acc] = {}
			if uniprot_resid not in uniprot_pdb_resids[uniprot_acc]:
				uniprot_pdb_resids[uniprot_acc][uniprot_resid] = set()
			uniprot_pdb_resids[uniprot_acc][uniprot_resid].add((chain, pdb_resid))
			if (chain, pdb_resid) not in pdb_uniprot_resids:
				pdb_uniprot_resids[(chain, pdb_resid)] = []
			if (uniprot_acc, uniprot_resid) not in pdb_uniprot_resids[(chain, pdb_resid)]:
				pdb_uniprot_resids[(chain, pdb_resid)].append((uniprot_acc, uniprot_resid))

#	print(indexed_pdb_uniprot_res_filename)
#	print("uniprot_pdb_resids", uniprot_pdb_resids)

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))
	print("\tDCA <<>> PDB\t", end='', flush=True)
	t = tf
	
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
			if new_convlist:
				new_mainlist.append(new_convlist)
		dca_pdb_resids[(pfam_acc, uniprot_acc)] = new_mainlist[:]

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))
	
	# List the residues with coordinates in each PDB chain
#	print("ALLOWED RESIDUES", pdbname)
	allowed_residues = {}
	for x in sorted(list(pdb_dca_resids.keys())):
		chain = x[0]
		resid = x[1]
		if chain not in allowed_residues:
			allowed_residues[chain] = set()
		allowed_residues[chain].add(resid)
#		print(pdbname, chain, resid)

#	print(pickle_filename)
	backmap_table = print_summary(pdb_dca_resids, pdb_uniprot_resids, backmap_filename)
	if len(target_pfam_accs) == 1:
		tcl_output_filename = results_folder + "/{0}_domain_visualization.tcl".format(target_pfam_accs[0][0])
	elif len(target_pfam_accs) == 2:
		tcl_output_filename = results_folder + "/{0}_{1}_domain_visualization.tcl".format(target_pfam_accs[0][0], target_pfam_accs[1][0])
	print("To visualize the domains: vmd -e {0}".format(tcl_output_filename))
	vmd_backmap_script(backmap_table, tcl_output_filename, pdb_path)
	pickle.dump((dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table), open(pickle_filename, 'wb'))

	print("")
	return dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table

def backmap_alignment(aln_path, pdb_path, chains, uniprot_id):
	"""Only the positions which correspond to DCA columns must be included in the alignment"""

	def gap_function(x, y):  # x is gap position in seq, y is gap length
		if y == 0:  # No gap
			return 0
		elif y == 1:  # Gap open penalty
			return -2
		return - (2 + y/4.0 + log(y)/2.0)


	detected_aln_format = None
	if aln_path.split('.')[-1] in ['fasta, fa']:
		detected_aln_format = 'fasta'
	elif aln_path.split('.')[-1] in ['stockholm, sto']:
		detected_aln_format = 'stockholm'
	else:
		with open(aln_path) as aln_f:
			for line in aln_f:
				if not line:
					continue
				if line.startswith('>'):
					detected_aln_format = 'fasta'
				elif line.startswith('# STOCKHOLM') or line.startswith('#=GF ID'):
					detected_aln_format = 'stockholm'
				break

	aln = []
	if not detected_aln_format:
		print('Error (align format not detected): could not detect format of alignment {0}'.format(aln_path))
	elif detected_aln_format == 'fasta':
		with open(aln_path) as aln_f:
			tmp_str = ''
			for line in aln_f:
				if not line:
					continue
				if line.startswith('>'):
					if tmp_str:
						aln.append(list(tmp_str))
				else:
					tmp_str += line.strip()
		if tmp_str:
			aln.append(list(tmp_str))
			
	elif detected_aln_format == 'stockholm':
		with open(aln_path) as aln_f:
			for line in aln_f:
				if not line:
					continue
				if not line.startswith('#') and line.strip() != '//':
					aln.append(list(line.split()[1]))

	carr = np.array(aln)
	del aln

	# Remove columns for which the first sequence have a gap
	del_col = []
	for i in carr.shape[1]:
		if carr[0,1].lower() == carr[0,1]:	# Removes small letters, '.' and '-': correct!
			del_col.append(i)
	
	np.delete(carr, del_col, axis=1)
	dca_model_length = corr.shape[1]

	str_seq = ''.join([x for x in carr[0,:]])	# All these characters must be mapped
	pdb_seqs = extract_sequence_from_pdb(pdb_path)

	for ch in pdb_seqs:
		str_pdb_aln = Bio.pairwise2.align.localdc(str_seq, pdb_seqs[ch], matlist.blosum62, gap_function, gap_function)



