from support import *
import matplotlib.pyplot as plt

def ppv(interaction_filename, dca_filename, dist_thr, distance_type, results_folder):
	distance_types = ['whole_res', 'sidechain', 'CA']

	if distance_type not in distance_types:
		print("ERROR: the chosen distance type is not valid: {0}".format(distance_type))
		exit(1)

	print("\nPPV:")
	print("distance type: {0}".format(distance_type))
	print("distance threshold: {0}\n".format(dist_thr))

	d = {}
	with open(interaction_filename) as interaction_file:
		for line in interaction_file:
			if not line:
				continue
			fields = line.split()
			if fields[5+distance_types.index(distance_type)] == 'None':
				d[(int(fields[0]), int(fields[1]))] = None
			else:
				d[(int(fields[0]), int(fields[1]))] = float(fields[6+distance_types.index(distance_type)])

	results_filename = results_folder + os.path.basename(dca_filename[:-4]) + 'ppv.txt'
	with open(dca_filename) as dca_file:
		n_contacts = 0
		ppv_score = []
		with open(results_filename, 'w') as results_file:
			for nl, line in enumerate(dca_file):
				if not line:
					continue
				fields = line.split()
				if (int(fields[0]), int(fields[1])) not in d or type(d[(int(fields[0]), int(fields[1]))]) == type(None):
					continue
				elif d[(int(fields[0]), int(fields[1]))] <= dist_thr:
					contact = True
					n_contacts += 1
				else:
					contact = False
				results_file.write("{0:6} {1:6} {2:10.3f} {3:10.3f} {4:7} {5:5} {6:10.3f}\n".format(fields[0], fields[1], fields[2], d[(int(fields[0]), int(fields[1]))], contact, n_contacts, n_contacts/(nl+1)))
				if nl < max(20,min(100,len(dca_file)/10)):
					print("{0:6} {1:6} {2:10.3f} {3:10.3f} {4:7} {5:5} {6:10.3f}".format(fields[0], fields[1], fields[2], d[(int(fields[0]), int(fields[1]))], contact, n_contacts, n_contacts/(nl+1)))
					ppv_score.append(n_contacts/(nl+1))

	results_figname = results_folder + os.path.basename(dca_filename[:-4]) + 'ppv.png'
	ppv_x = [x for x in range(1, len(ppv_score)+1)]
	plt.subplot(111)
	plt.plot(ppv_x, ppv_score)
	plt.savefig(results_figname)

	print("See figure relative to these data at {0}".figure(results_figname))
#	plt.show()
