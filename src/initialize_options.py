from support import *

def initialize_options(version, rel_path=False):
	if not version or version == 'None':
		print("ERROR: specify the version of the Pfam database used to create the DCA predictions (it's important!)")
		exit(1)

	print("\nPfam database version: {0}\n".format(version))

	options = {}
	if rel_path:
		this_path = "src/"
	else:
		this_path = os.path.dirname(os.path.abspath(__file__)) + '/'


	# Hardcoded paths
	options['src']                              = this_path
	options['support_files']                    = options['src']                    + 'support_files/'
	options['main']                             = options['src']                    + '../'
	options['db']                               = options['main']                   + 'db/'
	options['external_resources']               = options['db']                     + "external_resources/"
	options['pfam_version_main']                = options['external_resources']     + "Pfam_{0}/".format(version)
	options['database_files_relpath']           = options['pfam_version_main']      + "database_files/"            # Folder where the external files from the Pfam database are kept
	options['pdb_pfam_filename']                = options['database_files_relpath'] + "pdb_pfamA_reg.txt.gz"          # Register cennecting PDB and Pfam information
	options['pdb_uniprot_res_filename']         = options['database_files_relpath'] + "pdb_residue_data.txt.gz"       # Register with residue-by-residue details
	options['indexed_pdb_uniprot_res_folder']   = options['database_files_relpath'] + "pdb_residue_data/"
	options['indexed_pdb_uniprot_res_index']    = options['database_files_relpath'] + ".pdb_residue_data_index.txt"
	options['pfam_uniprot_filename']            = options['database_files_relpath'] + "uniprot_reg_full.txt.gz"       # Register connecting UniProt and Pfam information
	options['pfam_uniprot_stockholm_aln']       = options['database_files_relpath'] + "Pfam-A.full.uniprot.gz"
	options['pfam_pdbmap']                      = options['database_files_relpath'] + "pdbmap.gz"
	options['interaction_database']             = options['db']                     + "interaction_database/"
	options['pfam_pfam_filename']               = options['interaction_database']   + "PfamPfam_contacts_through_all_PDBs_table.txt"
	options['cache']                            = options['main']                   + ".cache/"
	options['pdb_files_ext_path']               = options['cache']                  + "PDB/"
	options['pfam_uniprot_stockholm_relpath']   = options['cache']                  + "alns_stockholm/"
	options['examples']                         = options['main']                   + 'examples/'
	options['results_folder']                   = options['examples']               + "results/"                   # Folder where results are put

	# Hardcoded options
	options['indexing']                       = True
	for p in [options['external_resources'], options['cache'], options['pdb_files_ext_path'], options['pfam_uniprot_stockholm_relpath'], options['results_folder']]:
		if not os.path.exists(p):
			os.mkdir(p)
	return options
