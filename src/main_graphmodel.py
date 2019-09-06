import matches
import backmap
import interactions
import mindistance
import networkx as nx
import matplotlib.pyplot as plt
from support import *


def search_pfam_for_uniprot(pfam_uniprot_stockholm_relpath, pfam_acc, ann_uniprot_acc):
	print(ann_uniprot_acc)
	msa_type = "uniprot"
	alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
	pfam_uniprot_stockholm_filename = pfam_uniprot_stockholm_relpath + pfam_acc + '_' + msa_type + '.stockholm'
	text = subprocess.run(['grep', '^{0}.'.format(ann_uniprot_acc[0]), pfam_uniprot_stockholm_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
	unpi, unpe = ann_uniprot_acc[1]
	for line in text:
		fields = line.split()
		i, e = [int(x) for x in fields[0].split('/')[1].split('-')]
		if i <= unpi and e >= unpe:
			return "".join([x for x in fields[1] if x in alphabet])


def calculate_hmm_seqID(seq1, seq2):
#	print(seq1)
#	print(seq2)
	if len(seq1) != len(seq2):
		print("ERROR: the two sequences do not have the same length")
		exit(1)
	if len(seq1)*len(seq2) == 0:
		print("ERROR: one of the two sequences is void")
		exit(1)

	n = 0
	hit = 0
	for i in range(len(seq1)):
		if seq1[i].upper() != seq1[i] or seq2[i].upper() != seq2[i]:
			continue
		if seq1[i] == seq2[i]:
			if seq1[i] == '-' or seq1[i] == '.':
				continue
			hit += 1
		n += 1
	return(hit/n)
	

def parallel_submission_routine(data):
	mind_pdbs, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, version, msa_type, force_download, pdb_files_ext_path, inpfam, inpfam1, inpfam2, inch1, inch2, pdb_pfam_filename, results_folder, cache_folder, self_inter = data
	main_backmap_table = {}
	interaction_filenames = set()
	for pdbname in mind_pdbs:
		print("\nPDB: ", pdbname, '-'*150)
		pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
		try:
			pfam_in_pdb = matches.calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename)
		except:
			failed_pdbs.add(pdbname)
			continue

		# Download PDB if needed
		pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
		if not os.path.exists(pdb_path):
			if not os.path.exists(pdb_files_ext_path):
				print("ERROR: PDB path not found")
				exit(1)
			download_pdb(pdbname, pdb_files_ext_path)
		if not os.path.exists(pdb_path):
			print("ERROR: PDB " + pdbname + " not found")
			exit(1)

		bundle = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
		if not bundle:
			failed_pdbs.add(pdbname)
			continue
		dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle
		main_backmap_table[pdbname] = backmap_table

		int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=self_inter)	
		interaction_filenames |= int_filenames
	return interaction_filenames, main_backmap_table


def calculate_connected_components(edge_labels, edge_weights):
	def rec(ccomp, neigh, v):
		for nv in neigh[v]:
			if nv not in ccomp:
				ccomp.append(nv)
				ccomp = rec(ccomp, neigh, nv)
		return ccomp


	vertices = []
	neighbors = {}
	for ie, edge_label in enumerate(edge_labels):
		for iv, v in enumerate(edge_label):
#			print(edge_label, iv, v)
			if v not in vertices:
				vertices.append(v)
				neighbors[v] = []
			if edge_weights[ie] == 0:
				continue
			if edge_label[1-iv] not in neighbors[v]:
				neighbors[v].append(edge_label[1-iv])
	visited_v = []
	conn_comps = []
	print("EDGES", [x for x in list(zip(edge_labels, edge_weights)) if x[1]>0])
	for v in vertices:
		print("VERTEX", v)
		if v not in visited_v:
			print("NOT VISITED")
			conn_comp = [v]
			visited_v.append(v)
			conn_comp = rec(conn_comp, neighbors, v)
			print("CONN_COMP", conn_comp)
			conn_comps.append(conn_comp)
			print("CONN_COMPS_TOT", conn_comps)
#			if len(conn_comp) > 1:
			visited_v += conn_comp#[1:]
			print("VISITED_V", visited_v)

	return conn_comps, vertices


def main_graphmodel(options):
	mindist = options['min_dist']
	inpfam = options['inpfam']
	inpfam1 = options['inpfam1']
	inpfam2 = options['inpfam2']
	pfam_pfam_filename = options['pfam_pfam_filename']
	pfam_pdbmap = options['pfam_pdbmap']
	pdb_pfam_filename = options['pdb_pfam_filename']
	pdb_uniprot_res_filename = options['pdb_uniprot_res_filename']
	indexed_pdb_uniprot_res_folder = options['indexed_pdb_uniprot_res_folder']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	pdb_files_ext_path = options['pdb_files_ext_path']
	results_folder = options['results_folder']
	cache_folder = options['cache']
	dca_filename = options['dca_filename']
	only_intra = options['only_intra']
	only_inter = options['only_inter']
	find_str = options['find_structures']
	pdb_uniprot_res_index_filename = options['indexed_pdb_uniprot_res_index']
	check_architecture = options['check_architecture']
	version = options['pfam_version']
	nprocs = options['nprocesses']
	inch1 = ''
	inch2 = ''

	dist1_filename = results_folder + "{0}_distances.txt".format(inpfam1)
	dist2_filename = results_folder + "{0}_distances.txt".format(inpfam2)
	distances1 = []
	distances1_labels = []
	distances2 = []
	distances2_labels = []
	if not (os.path.exists(dist1_filename) and os.path.exists(dist2_filename)):
		text = []
		pdbnames = []
		pdbtypes = [[], []]
		for i, inpfam in enumerate([inpfam1, inpfam2]):
			text.append(subprocess.run(["zgrep {0} {1} | awk '{{print substr($1,1,4)}}'".format(inpfam, pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n'))
			pdbnames.append(sorted(list(set([x.strip().lower() for x in text[i] if x.strip()])))) # DEBUG LEVA IL [:10] !!!!!!!!!!!!!!!!!!!!!!!!
			for pdbname in pdbnames[i]:
				textu = subprocess.run(["zgrep {0} {1} | awk '{{print substr($4,1,7)}}'".format(pdbname.upper(), pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
				textu = [x for x in textu if x.strip()]
				if not textu:
					print("Empty")
					exit(1)
				if len(textu) == 1:
					if textu[0] == inpfam:
						pdbtype = "single"
					else:
						print("How is this possible?")
						print(text)
						print(textu)
						exit(1)
				else:
					homogeneous = True
					has_partner = False
					for pfam_found in textu:
						if pfam_found != inpfam:
							homogeneous = False
						if pfam_found == [inpfam1, inpfam2][1-i]:
							has_partner = True
					if homogeneous:
						pdbtype = "homogeneous"
					elif not has_partner:
						pdbtype = "heterogeneous"
					else:
						pdbtype = "partner"
				pdbtypes[i].append(pdbtype)
		intersection_pdbs = sorted(list(set(pdbnames[0]) & set(pdbnames[1])))
	
		print("\n\nPDBs with both partners:")
		print("".join([x+" " for x in intersection_pdbs]))
		print("\n\nPDBS with Pfam {0}:".format(inpfam1))
		print("".join([x+" "+pdbtypes[0][i]+"\n" for i,x in enumerate(pdbnames[0])]))
		print("\n\nPDBS with Pfam {0}:".format(inpfam2))
		print("".join([x+" "+pdbtypes[1][i]+"\n" for i,x in enumerate(pdbnames[1])]))
	
		type0pdbs = [[x for i,x in enumerate(pdbnames[0]) if pdbtypes[0][i]==0], [x for i,x in enumerate(pdbnames[1]) if pdbtypes[1][i]==0]]
	
		with open("pdbnames_{0}.txt".format(inpfam1), "w") as f:
			for pdbname in sorted(list(set(pdbnames[0]) - set(intersection_pdbs))):
				f.write("{0}\n".format(pdbname))
		with open("pdbnames_{0}.txt".format(inpfam2), "w") as f:
			for pdbname in sorted(list(set(pdbnames[1]) - set(intersection_pdbs))):
				f.write("{0}\n".format(pdbname))
	
		failed_pdbs = set()	
		master_interactions = {}
		main_backmap_table = {}
		master_interactions_1 = {}
		main_backmap_table_1 = {}
		master_interactions_2 = {}
		main_backmap_table_2 = {}
		interaction_filenames = set()
		interaction_filenames_1 = set()
		interaction_filenames_2 = set()
		for pdbname in intersection_pdbs:
			print("\nPDB: ", pdbname, '-'*150)
			pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
			try:
				pfam_in_pdb = matches.calculate_matches(pdbname, '', inpfam1, inpfam2, pdb_pfam_filename)
			except:
				failed_pdbs.add(pdbname)
				continue
			print(failed_pdbs)
	
			# Download PDB if needed
			pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
			if not os.path.exists(pdb_path):
				if not os.path.exists(pdb_files_ext_path):
					print("ERROR: PDB path not found")
					exit(1)
				download_pdb(pdbname, pdb_files_ext_path)
			if not os.path.exists(pdb_path):
				print("ERROR: PDB " + pdbname + " not found")
				exit(1)

			bundle = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
			if not bundle:
				failed_pdbs.add(pdbname)
				continue
			dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle
			main_backmap_table[pdbname] = backmap_table
	
			int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=False)	
			interaction_filenames |= int_filenames
	
			for int_filename in int_filenames:
				with open(int_filename) as int_file:
					for line in int_file:
						fields = line.split()
						if fields[2] == 'None':
							dca1, dca2, ch1, resid1, ch2, resid2, dist = int(fields[0]), int(fields[1]), fields[2], fields[3], fields[4], fields[5], 1000000.0
						else:
							dca1, dca2, ch1, resid1, ch2, resid2, dist = int(fields[0]), int(fields[1]), fields[2], int(fields[3]), fields[4], int(fields[5]), float(fields[7])
						if (dca1, dca2) not in master_interactions:
							master_interactions[(dca1, dca2)] = []
						master_interactions[(dca1, dca2)].append((pdbname, (ch1, resid1), (ch2, resid2), dist))
	
	
		for pdbname in sorted(list(set(pdbnames[0]) - set(intersection_pdbs))):
			pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
			if not os.path.exists(pdb_path):
				if not os.path.exists(pdb_files_ext_path):
					print("ERROR: PDB path not found")
					exit(1)
				download_pdb(pdbname, pdb_files_ext_path)
			if not os.path.exists(pdb_path):
				print("ERROR: PDB " + pdbname + " not found")
				exit(1)
	
	
			pfam_in_pdb_1 = matches.calculate_matches(pdbname, inpfam1, '', '', pdb_pfam_filename)
	
			bundle = backmap.backmap_pfam(pfam_in_pdb_1, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
			if not bundle:
				failed_pdbs.add(pdbname)
				continue
			dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle
			main_backmap_table_1[pdbname] = backmap_table
	
			int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb_1, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=True)	
			interaction_filenames_1 |= int_filenames
	
			for int_filename in int_filenames:
				with open(int_filename) as int_file:
					for line in int_file:
						fields = line.split()
						if fields[2] == 'None':
							dca1, dca2, ch1, resid1, ch2, resid2, dist = int(fields[0]), int(fields[1]), fields[2], fields[3], fields[4], fields[5], 1000000.0
						else:
							dca1, dca2, ch1, resid1, ch2, resid2, dist = int(fields[0]), int(fields[1]), fields[2], int(fields[3]), fields[4], int(fields[5]), float(fields[7])
						if (dca1, dca2) not in master_interactions_1:
							master_interactions_1[(dca1, dca2)] = []
						master_interactions_1[(dca1, dca2)].append((pdbname, (ch1, resid1), (ch2, resid2), dist))
	
	
		offset = sorted(list(master_interactions.keys()))[-1][0]
		for pdbname in sorted(list(set(pdbnames[1]) - set(intersection_pdbs))):
			pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
			if not os.path.exists(pdb_path):
				if not os.path.exists(pdb_files_ext_path):
					print("ERROR: PDB path not found")
					exit(1)
				download_pdb(pdbname, pdb_files_ext_path)
			if not os.path.exists(pdb_path):
				print("ERROR: PDB " + pdbname + " not found")
				exit(1)
	
	
			pfam_in_pdb_2 = matches.calculate_matches(pdbname, inpfam2, '', '', pdb_pfam_filename)
	
			bundle = backmap.backmap_pfam(pfam_in_pdb_2, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
			if not bundle:
				failed_pdbs.add(pdbname)
				continue
			dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle
			main_backmap_table_2[pdbname] = backmap_table
	
	
			int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb_2, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=True)	
			interaction_filenames_2 |= int_filenames
	
			for int_filename in int_filenames:
				with open(int_filename) as int_file:
					for line in int_file:
						fields = line.split()
						if fields[2] == 'None':
							dca1, dca2, ch1, resid1, ch2, resid2, dist = int(fields[0]), int(fields[1]), fields[2], fields[3], fields[4], fields[5], 1000000.0
						else:
							dca1, dca2, ch1, resid1, ch2, resid2, dist = int(fields[0]), int(fields[1]), fields[2], int(fields[3]), fields[4], int(fields[5]), float(fields[7])
						if (dca1+offset, dca2+offset) not in master_interactions_2:
							master_interactions_2[(dca1+offset, dca2+offset)] = []
						master_interactions_2[(dca1+offset, dca2+offset)].append((pdbname, (ch1, resid1), (ch2, resid2), dist))
			
		interaction_dist_thr = 8.0
		interaction_list = []
		interaction_strength = []
		Ninter = len(master_interactions[list(master_interactions.keys())[0]])
		for dca1, dca2 in sorted(list(master_interactions.keys())):
			for i in range(Ninter):
				if master_interactions[(dca1, dca2)][i][3] <= interaction_dist_thr:
					if (dca1, dca2) not in interaction_list:
						interaction_list.append((dca1, dca2))
						interaction_strength.append(1/Ninter)
					else:
						interaction_strength[-1] += 1/Ninter
		
		print(interaction_list)
		print(interaction_strength)
	
		interaction_list_1 = []
		interaction_strength_1 = []
		Ninter = len(master_interactions_1[list(master_interactions_1.keys())[0]])
		vertices_1 = sorted(list(set([x[0] for x in interaction_list])))
		with open(results_folder + "{0}_distances.txt".format(inpfam1), "w") as f:
			for idca1, dca1 in enumerate(vertices_1):
				for idca2 in range(idca1+1, len(vertices_1)):
					dca2 = vertices_1[idca2]
					f.write("{0}\t{1}\t{2}\t{3}\t".format(dca1, dca2, interaction_strength[idca1], interaction_strength[idca2]))
					d1s = []
					for i in range(Ninter):
						d1s.append(master_interactions_1[(dca1, dca2)][i][3])
						f.write("{0}\t".format(master_interactions_1[(dca1, dca2)][i][3]))
					distances1_labels.append((dca1, dca2))
					distances1.append(d1s)
					f.write("\n")
	
	
		Ninter = len(master_interactions_2[list(master_interactions_2.keys())[0]])
		vertices_2 = sorted(list(set([x[1] for x in interaction_list])))
		with open(results_folder + "{0}_distances.txt".format(inpfam2), "w") as f:
			for idca1, dca1 in enumerate(vertices_2):
				for idca2 in range(idca1+1, len(vertices_2)):
					dca2 = vertices_2[idca2]
					f.write("{0}\t{1}\t{2}\t{3}\t".format(dca1, dca2, interaction_strength[idca1], interaction_strength[idca2]))
					d2s = []
					for i in range(Ninter):
						d2s.append(master_interactions_2[(dca1, dca2)][i][3])
						f.write("{0}\t".format(master_interactions_2[(dca1, dca2)][i][3]))
					distances2_labels.append((dca1, dca2))
					distances2.append(d2s)
					f.write("\n")

	else:
		vertices_1, vertices_2 = set(), set()
		with open(dist1_filename) as f:
			for line in f:
				if not line:
					continue
				fields = line.split()
				dca1, dca2, int_str1, int_str2 = int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3])
				vertices_1.add(dca1)
				vertices_1.add(dca2)
				distances1_labels.append((dca1, dca2))
				distances1.append([float(x) for x in fields[4:]])
		vertices_1 = sorted(list(vertices_1))
		with open(dist2_filename) as f:
			for line in f:
				if not line:
					continue
				fields = line.split()
				dca1, dca2, int_str1, int_str2 = int(fields[0]), int(fields[1]), float(fields[2]), float(fields[3])
				vertices_2.add(dca1)
				vertices_2.add(dca2)
				distances2_labels.append((dca1, dca2))
				distances2.append([float(x) for x in fields[4:]])
		vertices_2 = sorted(list(vertices_2))
	
	distances1 = np.array(distances1)
	distances2 = np.array(distances2)

	for cutoff_d in np.arange(2.0, 20.5, 0.5):
#		print(cutoff_d)
#		print(distances1 < cutoff_d)
		d1bool = ((distances1 < cutoff_d) & (distances1 >= 0))*1
		strengths1 = d1bool.sum(axis=1)
		print([(i,strengths1[i]) for i in range(len(strengths1)) if strengths1[i]>0])
		conn_comps1, vertices = calculate_connected_components(distances1_labels, strengths1)
		largest_size = sorted(conn_comps1, key= lambda x: -len(x))[0]
		print(inpfam1, cutoff_d, len(largest_size), largest_size)
#		print(inpfam1, sorted(conn_comps1, key= lambda x: -len(x)))

		d2bool = ((distances2 < cutoff_d) & (distances2 >= 0))*1
		strengths2 = d2bool.sum(axis=1)
		conn_comps2, vertices = calculate_connected_components(distances2_labels, strengths2)
		largest_size = sorted(conn_comps2, key= lambda x: -len(x))[0]
		print(inpfam2, cutoff_d, len(largest_size), largest_size)
#		print(inpfam2, sorted(conn_comps2, key= lambda x: -len(x)))


	# DCA PREDICTION
	# Find first pfam size
	offset = -1
	text = subprocess.run(["zcat {0}/{1}/{2}_uniprot_v{3}.stockholm.gz | tail -n10".format(pfam_uniprot_stockholm_relpath, inpfam1[:4], inpfam1, version)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
	for line in text:
		if not line or line.startswith("#") or line.startswith("//"):
			continue
		fields = line.split()
		offset = len([x for x in fields[1] if x == '-' or x.lower() != x])
		break
	if offset == -1:
		print("ERROR: no Pfam size found")
		exit(1)


	# read DCA
	recs = []
	with open(dca_filename) as dca_file:
		for line in dca_file:
			if not line or line.strip().startswith("#"):
				continue
			fields = line.split()
			recs.append((int(fields[0]), int(fields[1]), float(fields[2])))
#			print(dca_filename, int(fields[0]), int(fields[1]), float(fields[2]))
	recs = sorted(recs, key= lambda x: -x[2])

	# Retrieve edges	
	prediction_vertices_1 = []
	dca_distances1_labels, dca_distances2_labels, dca_distances1, dca_distances2 = [], [], [], []
	for dca1, dca2, score in recs:
		if dca1 <= offset and dca2 > offset:
			if dca1 not in prediction_vertices_1:
				prediction_vertices_1.append(dca1)
			dca_distances1_labels.append((dca1, dca2))
			if (dca1, dca2) in distances1_labels:
				dca_distances1.append(distances1[distances1_labels.index((dca1, dca2))])
			else:
				dca_distances1.append(1000000.0)
			if len(prediction_vertices_1) == 10:
				break
	dca_distances1 = np.array(distances1)
	dca_distances2 = np.array(distances2)

	for cutoff_d in np.arange(0.0, 20.5, 0.5):
		d1bool = ((distances1 < cutoff_d) & (distances1 >= 0))*1
		strengths1 = d1bool.sum(axis=1)
		conn_comps1, vertices = calculate_connected_components(dca_distances1_labels, strengths1)
		largest_size = sorted(conn_comps1, key= lambda x: -len(x))[0]
		print("DCA", inpfam1, cutoff_d, len(largest_size), largest_size)

		d2bool = ((distances2 < cutoff_d) & (distances2 >= 0))*1
		strengths2 = d2bool.sum(axis=1)
		conn_comps2, vertices = calculate_connected_components(dca_distances2_labels, strengths2)
		largest_size = sorted(conn_comps2, key= lambda x: -len(x))[0]
		print("DCA", inpfam2, cutoff_d, len(largest_size), largest_size)

	exit(1)	

	print("GRAPH A")	
	print(vertices_1)
	print(interaction_list_1)
	print(interaction_strength_1)

	number_of_first_vertices = 10
	vertices = sorted(list(set([x[0] for x in interaction_list]) | set([x[1] for x in interaction_list])))
	pair_list = sorted(list(zip(interaction_list_1, interaction_strength_1)), key= lambda x : x[1])
	drawn_vertices = []
	for x in pair_list:
		for e in x[0]:
			if e not in drawn_vertices and len(drawn_vertices) < 10:
				drawn_vertices.append(e)
		if len(drawn_vertices) == 10:
			break

	draw_interaction_list_1 = []
	draw_interaction_strength_1 = []
	draw_edges_1 = []
	for i, x in enumerate(interaction_list_1):
		if not(x[0] in drawn_vertices and x[1] in drawn_vertices):
			continue
		draw_interaction_list_1.append(x)
		draw_interaction_strength_1.append(interaction_strength_1[i])
		draw_edges_1.append((x[0], x[1], {'weight' : int(1+10*interaction_strength_1[i])}))
	G = nx.Graph()
	G.add_edges_from(draw_edges_1)
	plt.subplot(111)
	w = [G[u][v]['weight'] for u,v in G.edges]
	print("WIDTHS")
	print(w)
	nx.draw(G, with_labels=True, font_weight='bold', width=w)
	plt.savefig('prova.png')

	recs = []	
	with open(dca_filename) as dca_file:
		for line in dca_file:
			if not line or line.strip().startswith("#"):
				continue
			fields = line.split()
			recs.append((int(fields[0]), int(fields[1]), float(fields[2])))
#			print(dca_filename, int(fields[0]), int(fields[1]), float(fields[2]))
	recs = sorted(recs, key= lambda x: -x[2])


#	print("PREDICTION A:")
	prediction_vertices_1 = []
	for dca1, dca2, score in recs:
		if dca1 <= offset and dca2 > offset:
			if dca1 not in prediction_vertices_1:
				prediction_vertices_1.append(dca1)
			if len(prediction_vertices_1) == 10:
				break
	draw_interaction_list_1 = []
	draw_interaction_strength_1 = []
	draw_edges_1 = []
	for i, x in enumerate(interaction_list_1):
		if not(x[0] in prediction_vertices_1 and x[1] in prediction_vertices_1):
			continue
		draw_interaction_list_1.append(x)
		draw_interaction_strength_1.append(interaction_strength_1[i])
		draw_edges_1.append((x[0], x[1], {'weight' : int(1+10*interaction_strength_1[i])}))

	G = nx.Graph()
	G.add_edges_from(draw_edges_1)
	plt.clf()
	plt.subplot(111)
	w = [G[u][v]['weight'] for u,v in G.edges]
	print("WIDTHS")
	print(w)
	nx.draw(G, with_labels=True, font_weight='bold', width=w)
	plt.savefig('prova_pred.png')
	exit(1)


	interaction_list_2 = []
	interaction_strength_2 = []
	Ninter = len(master_interactions_2[list(master_interactions_2.keys())[0]])
	vertices_2 = sorted(list(set([x[1] for x in interaction_list])))
	for idca1, dca1 in enumerate(vertices_2):
		for dca2 in vertices_2[idca1+1:]:
			for i in range(Ninter):
				if master_interactions_2[(dca1, dca2)][i][3] <= interaction_dist_thr:
					if (dca1, dca2) not in interaction_list_2:
						interaction_list_2.append((dca1, dca2))
						interaction_strength_2.append(1/Ninter)
					else:
						interaction_strength_2[-1] += 1/Ninter

	print("GRAPH B")
	print(vertices_2)	
	print(interaction_list_2)
	print(interaction_strength_2)

	

	#for pdbname in 

	exit(1)

	if not dca_filename and not find_str:
		print("ERROR: mindist needs a precomputed plmdca filename")
	if find_str:
		mindist = 'all'
	
	mind_pdbs = []
	if mindist == 'all':
		if inpfam:
			self_inter = True
			text = subprocess.run(["zgrep {0} {1} | awk '{{print substr($1, 1, length($1)-1)}}'".format(inpfam, pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
			if check_architecture: 
#				print("grep {0} {1} | awk '$1==$2{{print $3}}'".format(inpfam, pfam_pfam_filename))
				textint = subprocess.run(["grep {0} {1} | awk '$1==$2{{print $3}}'".format(inpfam, pfam_pfam_filename)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
				textint = set(textint)
			for line in text:
				if not line:
					continue
				pdbname = line.strip().lower()
				if check_architecture:
					if pdbname not in textint:
						if pdbname not in mind_pdbs:
							mind_pdbs.append(pdbname)
					else:
#						print(pdbname, "removed")
						pass
				else:
					if pdbname not in mind_pdbs:
						mind_pdbs.append(pdbname)
		else:
			self_inter = False
			if check_architecture:
				text = subprocess.run(["grep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam1, pfam_pfam_filename, inpfam2)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
				pdbnames = sorted(list(set([x.strip().lower() for x in text if x.strip()])))
				new_pdbnames = []
				for pdbname in pdbnames:
					text1 = subprocess.run(["zgrep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam1, pfam_pdbmap, pdbname.upper())], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
					text2 = subprocess.run(["zgrep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam2, pfam_pdbmap, pdbname.upper())], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
					if text1.strip() and text2.strip() and pdbname not in new_pdbnames:
						mind_pdbs.append(pdbname)
			else:
				text1 = subprocess.run(["zgrep {0} {1} | awk '{{print substr($1,1,4)}}'".format(inpfam1, pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
				text2 = subprocess.run(["zgrep {0} {1} | awk '{{print substr($1,1,4)}}'".format(inpfam2, pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
				pdbnames1 = [x.strip().lower() for x in text1 if x.strip()]
				pdbnames2 = [x.strip().lower() for x in text2 if x.strip()]
				#print(pdbnames1, pdbnames2)
				mind_pdbs = sorted(list(set(pdbnames1) & set(pdbnames2)))
		mind_pdbs = sorted(mind_pdbs)
	else:
		if inpfam:
			self_inter = True
		else:
			self_inter = False
		if os.path.exists(mindist):
			with open(mindist) as mindist_file:
				for line in mindist_file:
					if not line:
						continue
					if line.strip() not in mind_pdbs:
						mind_pdbs.append(line.strip())
			mind_pdbs = sorted(mind_pdbs)
		else:
			print("ERROR: The file specifying the PDB names for the mindist analysis cannot be found")
			exit(1)

	structures_found_str = "".join([x+", " for x in mind_pdbs])[:-2]
	if find_str:
		print("\nPDBs found to contain the queried pfam(s):")
		print(structures_found_str+"\n")
		if inpfam:
			output_list_filename = results_folder + inpfam + "_structlist.txt"
		else:
			output_list_filename = results_folder + inpfam1 + "_" + inpfam2 + "_structlist.txt" 
		with open(output_list_filename, 'w') as output_list_file:
			for x in mind_pdbs:
				output_list_file.write(x+'\n')
		print("\nThis list was saved in {0}".format(output_list_filename))
		exit(1)
	else:
		print("\nConsidering the following PDBs:")
		print(structures_found_str+"\n")


	failed_pdbs = set()
	cluster = False
#	cluster = True
	if cluster and inpfam:
		sequences = {}
		maxi_backmap_table = []
		uniprots = {}
		for pdbname in mind_pdbs:#[:5]:	# DEBUG
			pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
			try:
				pfam_in_pdb = matches.calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename)
			except:
				failed_pdbs.add(pdbname)
				continue
	
			# Download PDB if needed
			pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
			if not os.path.exists(pdb_path) or os.stat(pdb_path).st_size == 0:
				if not os.path.exists(options['pdb_files_ext_path']):
					print("ERROR: PDB path not found")
					exit(1)
				download_pdb(pdbname, options['pdb_files_ext_path'])
			if not os.path.exists(pdb_path):
				print("ERROR: PDB " + pdbname + " not found")
				exit(1)

			bundle = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
			if not bundle:
				failed_pdbs.add(pdbname)
				continue
			dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle
			for line in backmap_table:
				unique_pfam_acc, pdb_corresp, unp_corresp = line
				pfam_acc = unique_pfam_acc.split('_')[0]
				if pfam_acc not in uniprots:
					uniprots[pfam_acc] = []
				if unp_corresp not in uniprots[pfam_acc]:
					uniprots[pfam_acc].append(unp_corresp)
				maxi_backmap_table.append([pdbname, unique_pfam_acc, pdb_corresp, unp_corresp])

		new_mind_pdbs = []
		for pdbname in mind_pdbs:
			if pdbname not in failed_pdbs:
				new_mind_pdbs.append(pdbname)
		mind_pdbs = new_mind_pdbs[:]

		print("HERE", uniprots)		


		seqs = {}
		pfam_set = {inpfam}     # MODIFY THIS FOR 2 PFAMS
		for pfam_acc in pfam_set:
			for unpi1, ann_uniprot_acc1 in enumerate(uniprots[pfam_acc]):
				seqs[(pfam_acc, ann_uniprot_acc1[0])] = search_pfam_for_uniprot(pfam_uniprot_stockholm_relpath, pfam_acc, ann_uniprot_acc1[0])
				print(pfam_acc, ann_uniprot_acc1[0], seqs[(pfam_acc, ann_uniprot_acc1[0])])

		for pfam_acc in pfam_set:
			hmm_seqID_mx = []
			for unpi1, ann_uniprot_acc1 in enumerate(uniprots[pfam_acc]):
				hmm_seqID_mx.append([])
				seq1 = seqs[(pfam_acc, ann_uniprot_acc1[0])]
				for unpi2, ann_uniprot_acc2 in enumerate(uniprots[pfam_acc]):
					if ann_uniprot_acc1 == ann_uniprot_acc2:
						hmm_seqID_mx[unpi1].append(1)
						continue
					seq2 = seqs[(pfam_acc, ann_uniprot_acc2[0])]
					hmm_seqID_mx[unpi1].append(calculate_hmm_seqID(seq1, seq2))

			N = len(uniprots[pfam_acc])
			hmm_seqIDs = np.ones((N, N))
			print(hmm_seqID_mx)
			for i in range(N):
				for j in range(N):
					print(i,j)
					hmm_seqIDs[i,j] = hmm_seqID_mx[i][j]

			hmm_seqIDs_labels = intrinsic_dimension_clustering(hmm_seqIDs)

#			dbscan(hmm_seqIDs)
#			exit(1)
		
#		for pdbname1 in mind_pdbs:
#			for pdbname2 in mind_pdbs:
#				for uniprot_acc1 in maxi_uniprot_pdb_resids[pdbname1]:
#					for uniprot_acc2 in maxi_uniprot_pdb_resids[pdbname2]:
							
			
	supermind_pdbs = []
	for i in range(nprocs):
		supermind_pdbs.append([])
	for i in range(len(mind_pdbs)):
		supermind_pdbs[i%nprocs].append(mind_pdbs[i])
	data = []
	for i in range(nprocs):
		data.append((supermind_pdbs[i], pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, version, msa_type, force_download, pdb_files_ext_path, inpfam, inpfam1, inpfam2, inch1, inch2, pdb_pfam_filename, results_folder, cache_folder, self_inter))

	pool = multiprocessing.Pool(processes=nprocs)
	pool_outputs = pool.map(parallel_submission_routine, data)
	pool.close()
	pool.join()

	main_backmap_table = {}
	interaction_filenames = set()
	for s in pool_outputs:
		interaction_filenames |= s[0]
		for x in s[1]:
			main_backmap_table[x] = s[1][x]

	"""
	main_backmap_table = {}
	interaction_filenames = set()
	for pdbname in mind_pdbs:
#		if cluster and pdbname in failed_pdbs:
#			continue
		print("\nPDB: ", pdbname, '-'*150)
		pdb_path = pdb_files_ext_path + pdbname.lower() + '.pdb'
		try:
			pfam_in_pdb = matches.calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename)
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

		dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, version, msa_type=msa_type, force_download=force_download)
		main_backmap_table[pdbname] = backmap_table

		int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=self_inter)	
		interaction_filenames |= int_filenames	
	"""

	new_mind_pdbs = []
	for x in mind_pdbs:
		if x in failed_pdbs:
			continue
		new_mind_pdbs.append(x)
	mind_pdbs = new_mind_pdbs

	interpolated_pdbs_string = "".join([x+", " for x in mind_pdbs])[:-2]
	print("\nPDBs effectively interpolated:")
	if not interpolated_pdbs_string:
		print("None. The algorithm will stop.")
		exit(1)
	else:
		print(interpolated_pdbs_string+"\n")

	if inpfam:
		mindistance.mindistance(mind_pdbs, inpfam, inpfam, only_intra, only_inter, results_folder, dca_filename, pfam_pdbmap, main_backmap_table, with_offset=False)
	else:
		mindistance.mindistance(mind_pdbs, inpfam1, inpfam2, only_intra, only_inter, results_folder, dca_filename, pfam_pdbmap, main_backmap_table)
