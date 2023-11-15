import sys
import os 
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../src/')

import initialize_options
import main_accessibility
import main_interactions
import main_mindistance
import main_backmap
import main_graphmodel
import interactions
from support import *

def main_parser():
	# Option parser
	parser = argparse.ArgumentParser()

	parser.add_argument('-pdb', '--pdbname', nargs='?')
	parser.add_argument('-a', '--accessibility', action='store_const', const='True', default='False')
	parser.add_argument('-adom', '--accessibilities_by_domain', action='store_const', const='True', default='False')
	parser.add_argument('-b', '--backmap', action='store_const', const='True', default='False')
	parser.add_argument('-arch', '--check_architecture', action='store_const', const='True', default='False')
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
	parser.add_argument('-distout', '--dist_filename', nargs='?')
	parser.add_argument('--only_distances', action='store_const', const='True', default='False')
	parser.add_argument('-v', '--pfam_version', nargs='?', default='32')
	parser.add_argument('-resolution', '--resolution_threshold', nargs='?', default='4.5')
	parser.add_argument('-restrict', '--restrict_comparison', nargs='?', default='')
	parser.add_argument('-out', '--output_path', nargs='?', default='None')
	parser.add_argument('-find_str', '--find_structures', action='store_const', const='True', default='False')
	parser.add_argument('-avg_dist', '--average_distance', action='store_const', const='True', default='False')
#	parser.add_argument('-no_indexing', '--no_indexing', action='store_const', const='True', default='False')
	parser.add_argument('--reset_cache', action='store_const', const='True', default='False')
	parser.add_argument('-no_cache', '--no_cache', action='store_const', const='True', default='False')
	parser.add_argument('-np', '--nprocesses', nargs='?')
	parser.add_argument('-graph', '--graphmodel', action='store_const', const='True', default='False')
	parser.add_argument('-compr', '--compress_distmx', action='store_const', const='True', default='False')

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
	parser.set_defaults(dist_filename = 'None')
	parser.set_defaults(nprocesses = '1')
	

	# Parse options
	parsed = parser.parse_args()
	options = initialize_options.initialize_options(version=parsed.pfam_version)
	for x in parsed.__dict__:
		options[x] = string_decode_element(parsed.__dict__[x], permissive=True, simple_version=True)

	if options['reset_cache']:
		subprocess.run("rm", "-f", "{0}/*".format(options['cache']))
		subprocess.run("rm", "-f", "{0}/.*".format(options['cache']))
	if options['dist_filename']:
		all_dist = interactions.compute_distances(options['pdbname'], options['dist_filename'], ch1=options['inch1'], ch2=options['inch2'])
		if options['only_distances']:
			exit(1)


	pfam_msa_types = {'uniprot' : 'uniprot', 'full' : 'full', 'rp75' : 'rp75'}
	if options['pfam_msa_type'] not in pfam_msa_types:
		print("ERROR: Pfam MSA type " + options['pfam_msa_type'] + " not defined")
		exit(1)
	options['msa_type'] = pfam_msa_types[options['pfam_msa_type']]

	if options['dist_filename']:
		all_dist = main_interactions.main_interactions(options)

	if options['output_path']:
		if os.path.exists(options['output_path']):
			if options['output_path'][-1] != '/':
				options['output_path'] += '/'
			options['results_folder'] = options['output_path'] + 'pfam2pdb_results_' + datetime.datetime.now().strftime("%Y%m%d") + "/"
			if not os.path.exists(options['results_folder']):
				os.mkdir(options['results_folder'])
		else:
			print("ERROR: the output path specified does not exist\n{0}".format(options['output_path']))

	if type(options['inpfam']) == type(None) and type(options['inpfam1']) == type(None) and type(options['inpfam2']) == type(None) and (not (options['backmap'])):
		print("ERROR: Specify at least one Pfam accession code")
		exit(1)


	printinfo = False
	inpfam_inpfam1_inpfam2 = []
	if options['inpfam'] == 'all':
		if not options['pdbname']:
			print("Error (incompatible options): -pf all needs -pdb")
			exit(1)
		if options['min_dist']:
			if options['min_dist'] != 'all':
				print("Error (incompatible options): -pf all can only support -mind all (with compulsory option -pdb)")
				exit(1)
			else:
				options['min_dist'] = '.automatic_min_dist.txt'
				with open(options['min_dist'], 'w') as f:
					f.write("{0}\n".format(options['pdbname']))
		text = subprocess.run(["zgrep {0} {1}".format(options['pdbname'].upper(), options['pfam_pdbmap'])], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')	
		set_of_pfams = set()
		for line in text:
			if not line:
				continue
			fields = line.split()
			for f in fields:
				if f[:2] == "PF":
					set_of_pfams.add(f[:-1])
					break
		for ippff, ppff1 in enumerate(sorted(list(set_of_pfams))):
			inpfam_inpfam1_inpfam2.append((None, ppff1, ppff1))	# same AND intra AND inter
			for ppff2 in sorted(list(set_of_pfams))[ippff+1:]:
				inpfam_inpfam1_inpfam2.append((None, ppff1, ppff2))
		printinfo = True
		print("Pfams present: {0}".format(set_of_pfams))
	else:
		for k in ['inpfam', 'inpfam1', 'inpfam2']:
			x = options[k]
			if type(x) != type(None) and (len(x) != 7 or x[:2] != "PF" or (not string_is_int(x[2:]))):
				print("ERROR: argument of option {0} is badly formatted".format(k))
				exit(1)
		inpfam_inpfam1_inpfam2.append((options['inpfam'], options['inpfam1'], options['inpfam2']))

	for i0, i1, i2 in inpfam_inpfam1_inpfam2:
		if printinfo:
			print("Structure: {0}".format(options['pdbname']))
			if i0:
				print("Pfam {0}".format(i0))
			else:
				print("Pfams {0}\t{1}".format(i1, i2))
		options['inpfam'], options['inpfam1'], options['inpfam2'] = i0, i1, i2
		if options['backmap'] == True:
			main_backmap.main_backmap(options)
		elif options['accessibility'] == True or options['accessibilities_by_domain']:
			main_accessibility.main_accessibility(options)
		elif type(options['min_dist']) != type(None) or options['find_structures'] == True:
			main_mindistance.main_mindistance(options)
		elif options['graphmodel'] == True:
			main_graphmodel.main_graphmodel(options)
		else:
			main_interactions.main_interactions(options)


if __name__ == "__main__":
	main_parser()
