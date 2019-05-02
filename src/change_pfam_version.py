import os
import sys
import initialize_options
from support import *


pfam_version = sys.argv[1]
if not string_is_int(pfam_version) or int(pfam_version) < 1:
	print("ERROR: version must be an integer > 0")
	exit(1)
pfam_version = pfam_version.strip()

options = initialize_options.initialize_options(pfam_version)

# Hardcoded paths
this_path = os.path.dirname(os.path.abspath(__file__)) + '/'
database_path = options['external_resources']
pfam_version_path = options['pfam_version_main']
database_filenames_filename = options['support_files'] + "database_filenames.txt"

if os.path.exists(pfam_version_path):
	print("Pfam version {0} is already present".format(pfam_version))
else:
	os.mkdir(options['pfam_version_main'])
	os.mkdir(options['database_files_relpath'])

	print("Downloading files from Pfam:")
	with open(database_filenames_filename) as database_filenames_file:
		for line in database_filenames_file:
			if line.startswith("#"):
				continue
			fields = line.split()
			database_filename = fields[0]
			if len(fields) == 2:
				database_folder = fields[1]
			else:
				database_folder = ""
			print("\t{0}".format(database_filename))
			os.system("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{0}.0/{2}/{1}.gz -O {3}/{1}.gz >/dev/null 2>&1".format(pfam_version, database_filename, database_folder, options['database_files_relpath']))

	os.system("{0}/decompression.sh {1} {2}".format(this_path, options['main'], pfam_version))

if options['indexing']:
	print("\nIndexing large files (this might take a while...\n")
	os.mkdir(options['indexed_pdb_uniprot_res_folder'])
	text = subprocess.run(["awk 'BEGIN{{fname=\"\"}}{{if (a[$1]!=1) {{if (fname!=\"\") {{close(fname)}}; x=substr($1,1,2); if (b[x]!=1) {{n[x]=1; b[x]=1}}; fname = \"{2}/\" x \"_{1}\"; print $1, fname, n[x]; a[$1]=1}}; print $0 >> fname; n[x]++}}' {0}".format(options['pdb_uniprot_res_filename'], os.path.basename(options['pdb_uniprot_res_filename']), options['indexed_pdb_uniprot_res_folder'])], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
	with open(options['indexed_pdb_uniprot_res_index'], 'w') as index_file:
		for line in text:
			index_file.write(line + '\n')

	with open(options['uniprot_stockholm_aln']) as aln_file:
		transcript = ""
		for line in aln_file:
			if not line:
				continue
			transcript += line + '\n'
			if line.startswith("#=GF AC"):
				fields = line.split()
				pfam_acc = fields[2].split(".")[0]
				current_filename = options[''] + current_filename + '_uniprot_v{0}.stockholm'.format(pfam_version)
			if line.startswith("//"):
				with open(current_filename, 'w') as out_file:
					out_file.write(transcript)
	os.remove(options['uniprot_stockholm_aln'])
