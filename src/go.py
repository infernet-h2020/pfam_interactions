from support import *

def strip_hydrogens(pdbname, pdb_path):

	conserved = []
	with open(pdb_path) as pdbf:
		for line in pdbf:
			if not line or not line.startswith("ATOM"):
				conserved.append(line)
				continue
			atname = line[11:16].strip()
			if atname[0] == "H":
				continue
			conserved.append(line)

	new_path = pdb_path
	if conserved:
		new_path = os.path.dirname(pdb_path) + '/' + 'nohydrogens_' + os.path.basename(pdb_path)
		print("WARNING! PDB FILE CONTAINS HYDROGENS. Creating file {0} which will be used for calculating interactions".format(new_path))
		with open (new_path, "w") as pdbf:
			for line in conserved:
				pdbf.write(line)
	return new_path


def compute_distances(pdbname, pdb_path, output_filename, distance_threshold, ch1='', ch2='', compress_distmx=False):
	distance = []
	sc_distance = []
	CA_distance = []
	new_tot_dist = {}
	pdb_path = strip_hydrogens(pdbname, pdb_path)
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbname, pdb_path)

	chain_pairs = {}
	with open(output_filename, 'w') as output_file:
		bb_names = {'N', 'O', 'OXT', 'C', 'CA'}
		chain_pairs = []
		for model in structure:
			for chain1 in model:
				ch1id = chain1.get_id()
				if ch1 and ch1 != ch1id:
					continue
				for chain2 in model:
					ch2id = chain2.get_id()
					if ch2 and ch2 != ch2id:
						continue
#					if (ch1id, ch2id) in chain_pairs or (ch2id, ch1id) in chain_pairs:
#						continue
					d = []
					for res1 in chain1:
						if res1.id[0].strip():
							continue
						for res2 in chain2:
							if res2.id[0].strip():
								continue
#							if ch1id == ch2id and res2.id[1] <= res1.id[1]:
#								continue

							if 'CA' in res1 and 'CA' in res2:
								CA_CA_d = np.linalg.norm(res1['CA'].get_vector() - res2['CA'].get_vector())
							else:
								CA_CA_d = -1

							CA_distance.append(CA_CA_d)
							if CA_CA_d > distance_threshold:
								sc_distance.append(CA_CA_d)
								distance.append(CA_CA_d)
								if compress_distmx:
									new_tot_dist[((ch1id, res1.id[1]), (ch2id, res2.id[1]))] = "MAX"
								else:
									output_file.write("{0}\t{1}\t{2}\t{3}\t{4:10.4f}\t{5:10.4f}\t{6:10.4f}\n".format(ch1id, res1.id[1], ch2id, res2.id[1], CA_CA_d, CA_CA_d, CA_CA_d))
									new_tot_dist[((ch1id, res1.id[1]), (ch2id, res2.id[1]))] = (distance[-1], sc_distance[-1], CA_distance[-1])
								continue

							mindist = 100000
							mindist_sc = 100000
							for a1 in res1:
								for a2 in res2:
									d = np.linalg.norm(a1.get_vector() - a2.get_vector())
									if d < mindist:
										mindist = d
									if d < mindist_sc and (a1.id not in bb_names) and (a2.id not in bb_names):
										mindist_sc = d
							if mindist == 100000:
								mindist = -1
							if mindist_sc == 100000:
								mindist_sc = -1
							distance.append(mindist)
							sc_distance.append(mindist_sc)
							new_tot_dist[((ch1id, res1.id[1]), (ch2id, res2.id[1]))] = (distance[-1], sc_distance[-1], CA_distance[-1])
							if not compress_distmx:
								output_file.write("{0}\t{1}\t{2}\t{3}\t{4:10.4f}\t{5:10.4f}\t{6:10.4f}\n".format(ch1id, res1.id[1], ch2id, res2.id[1], mindist, mindist_sc, CA_CA_d))
#			print(deep_getsizeof(distance, set())/(1024**2))
#			print(deep_getsizeof(sc_distance, set())/(1024**2))
#			print(deep_getsizeof(CA_distance, set())/(1024**2))
#			print(deep_getsizeof(new_tot_dist, set())/(1024**2))
#			print(deep_getsizeof(model, set())/(1024**2))
#			objgraph.show_refs(new_tot_dist, filename='sample-graph.png')
			break	# consider only first model

	if compress_distmx:
		os.remove(output_filename)
	return new_tot_dist


def compress_distance_matrix(dist_dict, distance_threshold):
#	print("Compress...")
	keys = sorted(list(set([k[0] for k in dist_dict])))
	linear_d = []
	set_list = []
	for ik1, k1 in enumerate(keys):
		for ik2, k2 in enumerate(keys):
			if type(dist_dict[(k1, k2)]) != str and dist_dict[(k1, k2)][2] <= distance_threshold and {k1, k2} not in set_list:
				linear_d.append(((k1, k2), dist_dict[(k1, k2)]))
				set_list.append({k1, k2})
	N = len(linear_d)
	np_lind = np.zeros(N, dtype=(float, 3))
	np_lind_idx = {}
	for i, x in enumerate(linear_d):
		a, b = x
		np_lind[i] = b
		idx1 = keys.index(a[0])
		idx2 = keys.index(a[1])
		np_lind_idx[(idx1, idx2)] = i
#	print("Done")
	return (np_lind, np_lind_idx, keys)


def decompress_distance_matrix(compressed_dmx):
#	print("Decompress...")
	np_lind, np_lind_idx, keys = compressed_dmx
	np_lind_idx_set = set(np_lind_idx)
	N = len(keys)
	dist_dict = {}
	for ik1, k1 in enumerate(keys):
		for ik2 in range(ik1, len(keys)):
			k2 = keys[ik2]
			if (ik1, ik2) in np_lind_idx_set:
				idx = np_lind_idx[(ik1, ik2)]
#				print(k1[0], k1[1], k2[0], k2[1], np_lind[idx])
				dist_dict[(k1, k2)] = np_lind[idx]
				dist_dict[(k2, k1)] = np_lind[idx]
#				np_lind_idx_set.discard((ik1, ik2))
#				np_lind_idx_set.discard((ik2, ik1))
			elif (ik2, ik1) in np_lind_idx_set:
				idx = np_lind_idx[(ik2, ik1)]
#				print(k1[0], k1[1], k2[0], k2[1], np_lind[idx])
				dist_dict[(k1, k2)] = np_lind[idx]
				dist_dict[(k2, k1)] = np_lind[idx]
#				np_lind_idx_set.discard((ik1, ik2))
#				np_lind_idx_set.discard((ik2, ik1))
			else:
#				print(k1[0], k1[1], k2[0], k2[1], "MAX")
				dist_dict[(k1, k2)] = "MAX"
				dist_dict[(k2, k1)] = "MAX"
#	print("Done")
	return dist_dict


def compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=False, allow_overwrite_dist=True, compress_distmx=False):
	print("Interactions:")
	t = time.time()
	bb_names = {'N', 'O', 'OXT', 'C'}
	distance_threshold = 20
	print("\tcalculating distances\t", end='', flush=True)
	pickle_filename = cache_folder + "." + pdbname + "_distances.pkl"
	pickle_compr_filename = cache_folder + "." + pdbname + "_compressed_distances.pkl"

	if ((not compress_distmx) and (not os.path.exists(pickle_filename))) or (compress_distmx and (not os.path.exists(pickle_compr_filename))):
		output_filename = results_folder + pdbname + "_distances.txt"
		tot_dist = compute_distances(pdbname, pdb_path, output_filename, distance_threshold, compress_distmx=compress_distmx)
		if compress_distmx:
			distmx_pkl = compress_distance_matrix(tot_dist, distance_threshold)
			pickle.dump(distmx_pkl, open(pickle_compr_filename, 'wb'))
		else:
			pickle.dump(tot_dist, open(pickle_filename, 'wb'))

	proto_residues_linear = []
	if ((not compress_distmx) and os.path.exists(pickle_filename)) or (compress_distmx and os.path.exists(pickle_compr_filename)):
		if compress_distmx:
			tot_dist = decompress_distance_matrix(pickle.load(open(pickle_compr_filename, 'rb')))
		else:
			tot_dist = pickle.load(open(pickle_filename, 'rb'))
		distance = []
		sc_distance = []
		CA_distance = []
		proto_residues_linear = []
		for x1, x2 in sorted(list(tot_dist.keys())):
			ch1, r1 = x1
			ch2, r2 = x2
			if ch1 in allowed_residues and r1 in allowed_residues[ch1] and ch2 in allowed_residues and r2 in allowed_residues[ch2]:
				if type(tot_dist[(x1, x2)]) == str and tot_dist[(x1, x2)] == "MAX":
					distance.append(9999)
					sc_distance.append(9999)
					CA_distance.append(9999)
				else:
					distance.append(tot_dist[(x1, x2)][0])
					sc_distance.append(tot_dist[(x1, x2)][1])
					CA_distance.append(tot_dist[(x1, x2)][2])
				if (ch1, r1) not in proto_residues_linear:
					proto_residues_linear.append((ch1, r1))

		residues_linear = proto_residues_linear
		print("(pickled!)\t", end='', flush=True)

	N = len(residues_linear)

	distance = np.array(distance).reshape((N, N))
	sc_distance = np.array(sc_distance).reshape((N, N))
	CA_distance = np.array(CA_distance).reshape((N, N))

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

							if val == 9999:
								interactions_file.write("{0:6}\t{1:6}\t{2:4}\t{3:4}\t{4:4}\t{5:4}\tMAX\tMAX\tMAX\t{6}\t{7}\n".format(i+1, j+1+offset, ac1, ar1, ac2, ar2, uniprot_resids1, uniprot_resids2))
							else:
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


if __name__ == "__main__":
	# pdbname, pdb_path, output_filename, distance_threshold
	compute_distances(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))
