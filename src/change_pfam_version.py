import os
import sys
from support import *

# Hardcoded paths
this_path = os.path.dirname(os.path.abspath(__file__)) + '/'
database_path = this_path + "../db/external_resources/"
database_filenames_filename = "support_files/database_filenames.txt"


pfam_version = sys.argv[1]
if not string_is_int(pfam_version) or int(pfam_version) < 1:
	print("ERROR: version must be an integer > 0")
	exit(1)
pfam_version = pfam_version.strip()
pfam_version_path = database_path + "Pfam_" + pfam_version + "/"	# WARNING: hardcoded name

if os.path.exists(pfam_version_path):
	print("Pfam version {0} is already present".format(pfam_version))
else:
	database_files_path = pfam_version_path + "database_files/"

	print(pfam_version_path)
	os.mkdir(pfam_version_path)
	os.mkdir(database_files_path)

	print("Downloading files from Pfam:")
	with open(database_filenames_filename) as database_filenames_file:
		for line in database_filenames_file:
			database_filename = line.strip()
			os.system("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{0}.0/database_files/{1}.gz -O {2}/{1}.gz".format(pfam_version, database_filename, database_files_path))
