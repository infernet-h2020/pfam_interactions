from support import *

def compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=False):
	print("Interactions:")
	t = time.time()
	bb_names = {'N', 'O', 'OXT', 'C'}
	distance_threshold = 30
	print("\tcalculating distances\t", end='', flush=True)
	pickle_filename = cache_folder + "." + pdbname + "_distances.pkl"

	proto_residues_linear = []
	unmapped_residues_linear = []
	if os.path.exists(pickle_filename):
		tot_dist, all_res = pickle.load(open(pickle_filename, 'rb'))
		distance = []
		sc_distance = []
		CA_distance = []
		proto_residues_linear = []
		for x1, x2 in sorted(list(tot_dist.keys())):
			ch1, r1 = x1
			ch2, r2 = x2
			if ch1 in allowed_residues and r1 in allowed_residues[ch1] and ch2 in allowed_residues and r2 in allowed_residues[ch2]:
				distance.append(tot_dist[(x1, x2)][0])
				sc_distance.append(tot_dist[(x1, x2)][1])
				CA_distance.append(tot_dist[(x1, x2)][2])
				if (ch1, r1) not in proto_residues_linear:
					proto_residues_linear.append((ch1, r1))
				if (ch2, r2) not in proto_residues_linear:
					proto_residues_linear.append((ch2, r2))
		proto_residues_linear = sorted(proto_residues_linear)

		allowed_residues_linear = []
		for ch in allowed_residues:
			for r in allowed_residues[ch]:
				allowed_residues_linear.append((ch, r))

		unmapped_residues_linear = sorted(list(set(allowed_residues_linear) - set(proto_residues_linear)))

#		N, residues_linear, distance, sc_distance, CA_distance = pickle.load(open(pickle_filename, 'rb'))
		print("(pickled!)\t", end='', flush=True)

	if not os.path.exists(pickle_filename) or unmapped_residues_linear:
		residues_linear = []
		distance = []
		sc_distance = []
		CA_distance = []
		parser = Bio.PDB.PDBParser(QUIET=True)
		structure = parser.get_structure(pdbname, pdb_path)

		new_tot_dist = {}
		for model in structure:
			for chain1 in allowed_residues:
				c1 = model[chain1]
				for r1 in c1:
					resid1 = r1.id[1]
					if r1.id[2].strip():
						continue
					if resid1 not in allowed_residues[chain1]:
						continue
					residues_linear.append((chain1, resid1))
					for chain2 in allowed_residues:
						c2 = model[chain2]
						if self_inter and chain2 != chain1:
							distance += [0]*len(allowed_residues[chain2])
							sc_distance += [0]*len(allowed_residues[chain2])
							CA_distance += [0]*len(allowed_residues[chain2])
							for resid2 in allowed_residues[chain2]:
								### ALT! POTREBBE ESSERCI UN BUG PIU' SERIO SE GLI ALLOWED RESIDUES SONO DIVERSI
								# IN PRATICA LE TRE RIGHE QUA SOPRA POTREBBERO PORTARE A UNA MATRICE SBAGLIATA (MA NON E' MAI SUCCESSO
								# PERCHE' C'ERA STO ALTRO BUG A MASCHERARE TUTTO!). INSOMMA, TOGLIERE! NON SI PUO' INSERIRE UNA RIGA DI
								# ZERI SENZA CONTROLLARE CHE CI SIANO NEL PDB! (MA FORSE CONTROLLO NEL BACKMAP?)
							for r2 in c2:
								resid2 = r2.id[1]
								if r2.id[2].strip():
									continue
								
							continue
						for r2 in c2:
							resid2 = r2.id[1]
							if r2.id[2].strip():
								continue
							if resid2 not in allowed_residues[chain2]:
								continue
							if (chain1 == chain2 and resid1 == resid2) or (inch1 and ({chain1, chain2} != {inch1, inch2})):
								distance.append(0)
								sc_distance.append(0)
								CA_distance.append(0)
								continue
							if (chain1, resid1) in proto_residues_linear and (chain2, resid2) in proto_residues_linear:
								distance.append(tot_dist[((chain1, resid1), (chain2, resid2))][0])
								sc_distance.append(tot_dist[((chain1, resid1), (chain2, resid2))][1])
								CA_distance.append(tot_dist[((chain1, resid1), (chain2, resid2))][2])
								continue
							CA_CA_d = np.linalg.norm(r1['CA'].get_vector() - r2['CA'].get_vector())
							CA_distance.append(CA_CA_d)
							if CA_CA_d > distance_threshold:
								sc_distance.append(CA_CA_d)
								distance.append(CA_CA_d)
								new_tot_dist[((chain1, resid1), (chain2, resid2))] = (CA_CA_d, CA_CA_d, CA_CA_d)
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
							new_tot_dist[((chain1, resid1), (chain2, resid2))] = (mindist, mindist_sc, CA_CA_d)
			break	# consider only first model
		N = len(residues_linear)

		distance = np.array(distance).reshape((N, N))
		sc_distance = np.array(sc_distance).reshape((N, N))
		CA_distance = np.array(CA_distance).reshape((N, N))

	if (not os.path.exists(pickle_filename)) or allow_overwrite_dist:
		if unmapped_residues_linear:
			for x in tot_dist:
				new_tot_dist[x] = tot_dist[x]
		pickle.dump((new_tot_dist, open(pickle_filename, 'wb'))

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
					if (not self_inter) and (not (pfam_acc1 == pfam_in_pdb[0][0] and pfam_acc2 == pfam_in_pdb[1][0])):	# Maintains the order of Pfams chosen at the beginning (in the command line!)
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
#					print(d, dSC, dCA)
					ausilia[(pfam_acc1, pfam_acc2)][(cmult1, cmult2)][dca_i1-1][dca_i2-1] = (c1, r1, c2, r2)

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))
	t = tf

#	print(ausilia)

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
			interactions_filename = results_folder + "{0}_{1}_{2}_{3}_{4}_{5}.txt".format(pfam_acc1, pfam_acc2, pdbname, cmult1, cmult2, imm)
			if pfam_acc1 == pfam_acc2 and cmult1 == cmult2:
				attr = 'same'
				offset = 0
			elif c1 == c2:
				attr = 'intra'
				offset = dca_model_length[pfam_acc1]
			else:
				attr = 'inter'
				offset = dca_model_length[pfam_acc1]
			if self_inter and attr != 'same':
				continue
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
	
							interactions_file.write("{0:6}\t{1:6}\t{2:4}\t{3:4}\t{4:4}\t{5:4}\t{6:12.6f}\t{7:12.6f}\t{8:12.6f}\t{9}\t{10}\n".format(i+1, j+1+offset, ac1, ar1, ac2, ar2, val, sc_val, CA_val, uniprot_resids1, uniprot_resids2))
			if (pfam_acc1, pfam_acc2) not in multiplicities:
				multiplicities[(pfam_acc1, pfam_acc2)] = []
			if (cmult1, cmult2) not in multiplicities[(pfam_acc1, pfam_acc2)]:
				multiplicities[(pfam_acc1, pfam_acc2)].append((cmult1, cmult2))
			multinum = multiplicities[(pfam_acc1, pfam_acc2)].index((cmult1, cmult2))
			interaction_summary.append((pfam_acc1, pfam_acc2, multinum, os.path.basename(interactions_filename), attr))
			interaction_filenames.add(interactions_filename)

	tf = time.time()
	print(time.strftime("%H:%M:%S", time.gmtime(tf-t)))

	for pfam_acc1, pfam_acc2, multinum, interactions_filename, attr in interaction_summary:
		print(pfam_acc1, pfam_acc2, multinum, interactions_filename, attr)

	if not interaction_filenames:
		print("")
	print("")
	return interaction_filenames
