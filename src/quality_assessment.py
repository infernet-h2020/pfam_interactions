import initialize_options
import main_ppv
from support import *

def main_parser():
	# Option parser
	parser = argparse.ArgumentParser()

	parser.add_argument('-int', '--interaction_filename', nargs='?', default='None')
	parser.add_argument('-dca', '--dca_filename', nargs='?', default='None')
	parser.add_argument('-thr', '--dist_thr', nargs='?', default='8')
	parser.add_argument('-dtype', '--distance_type', nargs='?', default='whole_res')

	# Parse options
	parsed = parser.parse_args()
	options = initialize_options.initialize_options(version='32')
#	print(parsed)
	for x in parsed.__dict__:
		options[x] = string_decode_element(parsed.__dict__[x], permissive=True, simple_version=True)
#	print(options)

	main_ppv.ppv(options['interaction_filename'], options['dca_filename'], options['dist_thr'], options['distance_type'])

if __name__ == "__main__":
	main_parser()
