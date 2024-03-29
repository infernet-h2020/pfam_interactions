import matches
import backmap
import interactions
from support import *

def main_interactions(options):
	pdbname = options['pdbname']
	inpfam = options['inpfam']
	inpfam1 = options['inpfam1']
	inpfam2 = options['inpfam2']
	pdb_pfam_filename = options['pdb_pfam_filename']
	pdb_uniprot_res_filename = options['pdb_uniprot_res_filename']
	indexed_pdb_uniprot_res_folder = options['indexed_pdb_uniprot_res_folder']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	inch1 = options['inch1']
	inch2 = options['inch2']
	results_folder = options['results_folder']
	cache_folder = options['cache']
	pdb_uniprot_res_index_filename = options['indexed_pdb_uniprot_res_index']
	version = options['pfam_version']
	dist_filename = options['dist_filename']
	only_distances = options['only_distances']
	compress_distmx = options['compress_distmx']
	resolution_threshold = options['resolution_threshold']

	# Download PDB if needed
	pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
	if not os.path.exists(pdb_path):
		if not os.path.exists(options['pdb_files_ext_path']) or os.stat(options['pdb_files_ext_path']).st_size == 0:
			print("ERROR: PDB path not found")
			print(options['pdb_files_ext_path'])
			exit(1)
		download_pdb(pdbname, options['pdb_files_ext_path'])
	if not os.path.exists(pdb_path):
		print("ERROR: PDB " + pdbname + " not found")
		exit(1)

	if not check_pdb_quality(pdb_path, resolution_threshold):
		exit(1)

	options['pdb_path'] = pdb_path

	if options['dist_filename']:
		all_dist = interactions.compute_distances(pdbname, pdb_path, dist_filename, ch1=inch1, ch2=inch2)
		if only_distances:
			exit(1)

	pfam_in_pdb = matches.calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename)

	bundle = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
	if not bundle:
		exit(1)
	dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle

	if inpfam and (not inpfam1) and (not inpfam2):
		self_inter = True
	else:
		self_inter = False
	interactions.compute_interactions(pdbname, pdb_path, pfam_in_pdb, pdb_uniprot_resids, uniprot_restypes, pdb_dca_resids, dca_model_length, allowed_residues, inch1, inch2, results_folder, cache_folder, self_inter=self_inter, compress_distmx=compress_distmx)
