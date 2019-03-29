import matches
import backmap
import interactions
import mindistance
from support import *

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
	check_architecture = True
	inch1 = ''
	inch2 = ''

	if not dca_filename and not find_str:
		print("ERROR: mindist needs a precomputed plmdca filename")
	if find_str:
		mindist = 'all'
	
	mind_pdbs = []
	if mindist == 'all':
		if inpfam:
			self_inter = True
			text = subprocess.run(["grep {0} {1} | awk '{{print substr($1, 1, length($1)-1)}}'".format(inpfam, pfam_pdbmap)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
			if check_architecture:
				print("grep {0} {1} | awk '$1==$2{{print $3}}'".format(inpfam, pfam_pfam_filename))
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
						print(pdbname, "removed")
		else:
			self_inter = False
			text = subprocess.run(["grep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam1, pfam_pfam_filename, inpfam2)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
			for line in text:
				if not line:
					continue
				pdbname = line.strip().lower()
				text1 = subprocess.run(["grep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam1, pfam_pdbmap, pdbname.upper())], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
				text2 = subprocess.run(["grep {0} {1} | grep {2} | awk '{{print $3}}'".format(inpfam2, pfam_pdbmap, pdbname.upper())], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
				if text1.strip() and text2.strip():
					if pdbname not in mind_pdbs:
						mind_pdbs.append(pdbname)
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

	### METTERLO QUI E' UN CASINO, DEVI FARLO 
#	if cluster:
#		if inpfam:
#			inpfams = [inpfams]
#		else:
#			inpfams = [inpfam1, inpfam2]
#		seqID_d = {}
#		seqID = np.ones((len(mind_pdbs), len(mind_pdbs)))
#		for i, pdbname1 in enumerate(mind_pdbs):
#			for j, pdbname2 in enumerate(mind_pdbs[i+1:]):
#				seqs_pdb1 = get_Pfam_seq(inpfams, pdbname1)
#				seqs_pdb2 = get_seq(pdbname2)
#				matrix = matlist.blosum62
#				alns = pairwise2.align.globalds(seq1, seq2, matrix, -2, -0.25)
#				seqID[i][j] = min(calculate_asymm_seqid(alns[0]))	# Take the minimum between the two seqIDs relative to the first alignment
#				seqID[j][i] = seqID[i][j]


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
		pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
		if not os.path.exists(pdb_path):
			if not os.path.exists(options['pdb_files_ext_path']):
				print("ERROR: PDB path not found")
				exit(1)
			download_pdb(pdbname, options['pdb_files_ext_path'])
		if not os.path.exists(pdb_path):
			print("ERROR: PDB " + pdbname + " not found")
			exit(1)

		print(pfam_in_pdb)

		dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pfam_uniprot_stockholm_relpath, cache_folder, msa_type=msa_type, force_download=force_download)

		int_filenames = interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=self_inter)	
		interaction_filenames |= int_filenames	

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
		mindistance.mindistance(mind_pdbs, inpfam, inpfam, only_intra, only_inter, results_folder, dca_filename, pfam_pfam_filename)
	else:
		mindistance.mindistance(mind_pdbs, inpfam1, inpfam2, only_intra, only_inter, results_folder, dca_filename, pfam_pfam_filename)
