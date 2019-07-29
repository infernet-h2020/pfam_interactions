import matches
import backmap
import interactions
import mindistance
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
	mind_pdbs, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, version, msa_type, force_download, pdb_files_ext_path, inpfam, inpfam1, inpfam2, inch1, inch2, pdb_pfam_filename, results_folder, cache_folder, self_inter, compress_distmx = data
	main_backmap_table = {}
	interaction_filenames = set()
	failed_pdbs = set()
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

		int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=self_inter, compress_distmx=compress_distmx)
		interaction_filenames |= int_filenames
	return interaction_filenames, main_backmap_table


def main_mindistance(options):
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
	compress_distmx = options['compress_distmx']
	inch1 = ''
	inch2 = ''

	if not dca_filename and not find_str:
		print("WARNING: mindist will not be ordered by a DCA score (DCA file is missing)")
	if find_str:
		mindist = 'all'
	
	mind_pdbs = []
	if mindist == 'all':
		if inpfam:
			self_inter = True
			text = subprocess.run(["zgrep {0} {1} | awk '{{print substr($1, 1, length($1)-1)}}'".format(inpfam, pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
#			print("zgrep {0} {1} | awk '{{print substr($1, 1, length($1)-1)}}'".format(inpfam, pfam_pdbmap))
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
#				print(pdbnames1, pdbnames2)
				mind_pdbs = sorted(list(set(pdbnames1) & set(pdbnames2)))
#				print(mind_pdbs)
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

			dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
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
		data.append((supermind_pdbs[i], pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, version, msa_type, force_download, pdb_files_ext_path, inpfam, inpfam1, inpfam2, inch1, inch2, pdb_pfam_filename, results_folder, cache_folder, self_inter, compress_distmx))

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
