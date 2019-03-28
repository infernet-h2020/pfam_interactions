#!/bin/python3

import initialize_options
import main_accessibility
import main_interactions
import main_mindistance
from support import *

def main_parser():
	# Option parser
	parser = argparse.ArgumentParser()

	parser.add_argument('-pdb', '--pdbname', nargs='?')
	parser.add_argument('-a', '--accessibility', action='store_const', const='True', default='False')
	parser.add_argument('-pf', '--inpfam', nargs='?')
	parser.add_argument('-pf1', '--inpfam1', nargs='?')
	parser.add_argument('-pf2', '--inpfam2', nargs='?')
	parser.add_argument('-c1', '--inch1', nargs='?')
	parser.add_argument('-c2', '--inch2', nargs='?')
	parser.add_argument('-ptype', '--pfam_msa_type', nargs='?')
	parser.add_argument('-vmd', '--vmd_path', nargs='?')
	parser.add_argument('-fd', '--force_download', nargs='?')
	parser.add_argument('-mind', '--min_dist', nargs='?')
	parser.add_argument('-dca', '--dca_filename', nargs='?')
	parser.add_argument('-intra', '--only_intra', action='store_const', const='True', default='False')
	parser.add_argument('-inter', '--only_inter', action='store_const', const='True', default='False')
	parser.add_argument('-v', '--pfam_version', nargs='?', default='32')
	parser.add_argument('-out', '--output_path', nargs='?', default='None')
	parser.add_argument('-find_str', '--find_structures', action='store_const', const='True', default='False')

	# Default values for optional arguments
	parser.set_defaults(pdbname = 'None')
	parser.set_defaults(accessibility = 'False')
	parser.set_defaults(inpfam = 'None')
	parser.set_defaults(inpfam1 = 'None')
	parser.set_defaults(inpfam2 = 'None')
	parser.set_defaults(inch1 = 'None')
	parser.set_defaults(inch2 = 'None')
	parser.set_defaults(pfam_msa_type = 'uniprot')
	parser.set_defaults(vmd_path = '/Applications/VMD\ 1.9.3.app/Contents/vmd/vmd_MACOSXX86')
	parser.set_defaults(force_download = 'False')
	parser.set_defaults(min_dist = 'None')
	parser.set_defaults(dca_filename = 'None')
	

	# Parse options
	parsed = parser.parse_args()
	options = initialize_options.initialize_options(version=parsed.pfam_version)
#	print(parsed)
	for x in parsed.__dict__:
		options[x] = string_decode_element(parsed.__dict__[x], permissive=True, simple_version=True)
#	print(options)


	pfam_msa_types = {'uniprot' : 'uniprot', 'full' : 'full', 'rp75' : 'rp75'}
	if options['pfam_msa_type'] not in pfam_msa_types:
		print("ERROR: Pfam MSA type " + options['pfam_msa_type'] + " not defined")
		exit(1)
	options['msa_type'] = pfam_msa_types[options['pfam_msa_type']]

	if options['output_path']:
		if os.path.exists(options['output_path']):
			options['results_folder'] = options['output_path'] + 'dca2pdb_results_' + datetime.datetime.now().strftime("%Y%m%d") + "/"
			if not os.path.exists(options['results_folder']):
				os.mkdir(options['results_folder'])
		else:
			print("ERROR: the output path specified does not exist\n{0}".format(options['output_path']))

	if type(options['inpfam']) == type(None) and type(options['inpfam1']) == type(None) and type(options['inpfam2']) == type(None):
		print("ERROR: Specify at least one Pfam accession code")
		exit(1)

	for k in ['inpfam', 'inpfam1', 'inpfam2']:
		x = options[k]
		if type(x) != type(None) and (len(x) != 7 or x[:2] != "PF" or (not string_is_int(x[2:]))):
			print("ERROR: argument of option {0} is badly formatted".format(k))
			exit(1)

	if options['accessibility'] == True:
		main_accessibility.main_accessibility(options)
	elif type(options['min_dist']) != type(None) or options['find_structures'] == True:
		main_mindistance.main_mindistance(options)
	else:
		main_interactions.main_interactions(options)


if __name__ == "__main__":
	main_parser()
