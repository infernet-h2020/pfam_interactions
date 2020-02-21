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
	indexed_pdb_uniprot_res_folder = options['indexed_pdb_uniprot_res_folder']
	pfam_uniprot_stockholm_relpath = options['pfam_uniprot_stockholm_relpath']
	msa_type = options['msa_type']
	force_download = options['force_download']
	vmd_path = options['vmd_path']
	results_folder = options['results_folder']
	cache_folder = options['cache']
	pdb_uniprot_res_index_filename = options['indexed_pdb_uniprot_res_index']
	version = options['pfam_version']
	accessibilities_by_domain = options['accessibilities_by_domain']
	resolution_threshold = options['resolution_threshold']
	
	if (not vmd_path):
		print("ERROR: SASA calculation needs VMD to run. Specify a path for VMD")
		exit(1)

	# Download PDB if needed
	pdb_path = options['pdb_files_ext_path'] + pdbname.lower() + '.pdb'
	if not os.path.exists(pdb_path):
		if not os.path.exists(options['pdb_files_ext_path']) or os.stat(options['pdb_files_ext_path']).st_size == 0:
			print("ERROR: PDB path not found")
			exit(1)
		download_pdb(pdbname, options['pdb_files_ext_path'])
	if not os.path.exists(pdb_path):
		print("ERROR: PDB " + pdbname + " not found")
		exit(1)

	if not check_pdb_quality(pdb_path, resolution_threshold):
		exit(1)

	options['pdb_path'] = pdb_path

	pfam_in_pdb = matches.calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename)

	bundle = backmap.backmap_pfam(pfam_in_pdb, pdbname, pdb_path, pdb_pfam_filename, pdb_uniprot_res_filename, indexed_pdb_uniprot_res_folder, pdb_uniprot_res_index_filename, pfam_uniprot_stockholm_relpath, cache_folder, results_folder, version, msa_type=msa_type, force_download=force_download)
	if not bundle:
		exit(1)
	dca_model_length, uniprot_restypes, uniprot_pdb_resids, pdb_uniprot_resids, dca_pdb_resids, pdb_dca_resids, allowed_residues, backmap_table = bundle

	SASA = accessibility.compute_accessibility(pdbname, pdb_path, pdb_dca_resids, backmap_table, allowed_residues, vmd_path, results_folder, accessibilities_by_domain)

	if inpfam1 and inpfam2:
		n1, n2 = 0, 0
		avg_SASA1 = -np.ones(dca_model_length[inpfam1])
		avg_SASA2 = -np.ones(dca_model_length[inpfam2])
		avg_SASA = -np.ones((dca_model_length[inpfam1], dca_model_length[inpfam2]))
		solvation_map_fig = results_folder + '/{0}_{1}_{2}_avg_solvent_accessibility.png'.format(pdbname, inpfam1, inpfam2)
		for pfam_acc in SASA:
			if inpfam1 in pfam_acc:
				n1 += 1
				for k in SASA[pfam_acc]:
					dca_i, _, _ = k
					if avg_SASA1[dca_i-1] == -1:
						avg_SASA1[dca_i-1] = 0
					avg_SASA1[dca_i-1] += SASA[pfam_acc][k][1]
			if inpfam2 in pfam_acc:
				n2 += 1
				for k in SASA[pfam_acc]:
					dca_i, _, _ = k
					if avg_SASA2[dca_i-1] == -1:
						avg_SASA2[dca_i-1] = 0
					avg_SASA2[dca_i-1] += SASA[pfam_acc][k][1]
		for i in range(avg_SASA1.shape[0]):
			for j in range(avg_SASA2.shape[0]):
				if avg_SASA1[i] == -1 or avg_SASA2[j] == -1:
					avg_SASA[i,j] = -0.5
				else:
					avg_SASA[i,j] = (avg_SASA1[i] * avg_SASA2[j]) / (n1 * n2)

		plt.imshow(avg_SASA, cmap='BrBG', interpolation='none', vmin=-1, vmax=1)
		plt.savefig(solvation_map_fig)
