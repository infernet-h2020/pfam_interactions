import backmap
from support import *

def main_backmap(options):
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
	vmd_path = options['vmd_path']
	results_folder = options['results_folder']
	cache_folder = options['cache']
	pdb_uniprot_res_index_filename = options['indexed_pdb_uniprot_res_index']
	version = options['pfam_version']
	complete_backmap = False
	resolution_threshold = options['resolution_threshold']
	
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

	if not check_pdb_quality(pdb_path, resolution_threshold):
		exit(1)

	options['pdb_path'] = pdb_path

	if inpfam:
		pfam_in_pdb = [(inpfam, "XXX")]
	elif (inpfam1 and inpfam2):
		pfam_in_pdb = [(inpfam1, "XXX"), (inpfam2, "XXX")]
	else:
		pfam_in_pdb = []
		complete_backmap = True


	bundle = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download, complete_backmap=complete_backmap)
	if not bundle:
		exit(1)
