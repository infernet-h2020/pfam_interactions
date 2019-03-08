import matches
import backmap
import accessibility
from support import *

def main_accessibility(options):
	pdbname = options['pdbname']
	inpfam = options['inpfam']
	inpfam1 = options['inpfam1']
	inpfam2 = options['inpfam2']
	pdb_pfam_filename = options['pdb_pfam_filename']
	pdb_uniprot_res_filename = options['pdb_uniprot_res_filename']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	vmd_path = options['vmd_path']
	results_folder = options['results_folder']
	cache_folder = options['cache']
	
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

	pfam_in_pdb = matches.calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename)

	dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_pfam_filename, pdb_uniprot_res_filename, pfam_uniprot_stockholm_relpath, cache_folder, msa_type=msa_type, force_download=force_download)

	SASA = accessibility.compute_accessibility(pdbname, pdb_path, pdb_dca_resids, allowed_residues, vmd_path, results_folder)
