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
	totspace = 0
	downs = []
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
			URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{0}.0/{2}/{1}.gz".format(pfam_version, database_filename, database_folder)
			destination = "{1}/{0}.gz".format(database_filename, options['database_files_relpath'])
			text = subprocess.run(["curl -sI {0} | grep -i Content-Length | awk '{{print $2}}' | {1}/byte_to_human.sh".format(URL, this_path)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
			totspace += int(subprocess.run(["curl -sI {0} | grep -i Content-Length".format(URL)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split("\n")[0].split()[1])
			num, unit = text.split('\n')[0].split()
			print("\t{0}\t{1} {2}".format(database_filename, num, unit))
			downs.append((URL, destination))
	totnum, totunit = subprocess.run(["echo {1} | {0}/byte_to_human.sh".format(this_path, totspace)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split("\n")[0].split()
	print("Total disk space required:\t{0} {1}".format(totnum, totunit))
	ans = query_yes_no("Proceed with download?", default="yes")
	if ans:
		print("Downloading (this may take a while...)")
		for URL, destination in downs:
			subprocess.run(["wget", URL, "-O", destination], stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))
	else:
		print("Installation interrupted")
		os.rmdir(options['database_files_relpath'])
		os.rmdir(options['pfam_version_main'])
		exit(1)

#	subprocess.run(["{0}/decompression.sh".format(this_path), options['main'], pfam_version], stdout=open("/dev/null", "w"), stderr=open("/dev/null", "w"))

if options['indexing']:
	print("\nIndexing large files (this might take a while...)\n")
	if not os.path.exists(options['indexed_pdb_uniprot_res_folder']):
		os.mkdir(options['indexed_pdb_uniprot_res_folder'])

	text = subprocess.run(["gzcat {0} | awk 'BEGIN{{fname=\"\"}}{{if (a[$1]!=1) {{if (fname!=\"\") {{close(fname)}}; x=substr($1,1,2); if (b[x]!=1) {{n[x]=1; b[x]=1}}; fname = \"{2}/\" x \"_{1}\"; print $1, fname, n[x]; a[$1]=1}}; print $0 >> fname; n[x]++}}'".format(options['pdb_uniprot_res_filename'], os.path.basename(options['pdb_uniprot_res_filename'])[:-3], options['indexed_pdb_uniprot_res_folder'])], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')

	with open(options['indexed_pdb_uniprot_res_index'], 'w') as index_file:
		for line in text:
			index_file.write(line + '\n')


	transcript = []
	for line in codecs.getreader('utf-8')(gzip.open(options['pfam_uniprot_stockholm_aln']), errors='replace'):
#	with gzip.open(options['pfam_uniprot_stockholm_aln'], 'rt', encoding='utf-8') as aln_file:
#		transcript = []
#		while True:
#			line = aln_file.readline()
#			if line:
		transcript.append(line)
		if line.startswith("#=GF AC"):
			fields = line.split()
			pfam_acc = fields[2].split(".")[0]
			current_filename = options['pfam_uniprot_stockholm_relpath'] + pfam_acc + '_uniprot_v{0}.stockholm'.format(pfam_version)

		if line.startswith("//"):
			print(current_filename)
			index_folder = options['pfam_uniprot_stockholm_relpath'] + pfam_acc[:4] + '/'
			if not os.path.exists(index_folder):
				os.mkdir(index_folder)
			current_filename = index_folder + pfam_acc + '_uniprot_v{0}.stockholm'.format(pfam_version)
			zipped_filename = current_filename + ".gz"
			if os.path.exists(zipped_filename):
				transcript = []
				continue
			with open(current_filename, 'w') as out_file:
				bl = 10000
				for b in range(int(len(transcript)/bl)+1):
					out_file.write("".join(transcript[b*bl:(b+1)*bl]))
#			print("write", zipped_filename)
			subprocess.run(["gzip", current_filename])
			transcript = []
#			else:
#				break
	os.remove(options['pfam_uniprot_stockholm_aln'])
