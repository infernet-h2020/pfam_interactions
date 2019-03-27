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
			database_filename = line.strip()
			print("\t{0}".format(database_filename))
			os.system("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{0}.0/database_files/{1}.gz -O {2}/{1}.gz 2>&1 > /dev/null".format(pfam_version, database_filename, options['database_files_relpath']))
