from support import *
import matplotlib.pyplot as plt

def ppv(interaction_filename, dca_filename, dist_thr, distance_type):
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

	with open(dca_filename) as dca_file:
		n_contacts = 0
		ppv_score = []
		for nl, line in enumerate(dca_file):
			if not line:
				continue
			fields = line.split()
			if type(d[(int(fields[0]), int(fields[1]))]) == type(None):
				continue
			elif d[(int(fields[0]), int(fields[1]))] <= dist_thr:
				contact = True
				n_contacts += 1
			else:
				contact = False
			print(fields[0], fields[1], fields[2], d[(int(fields[0]), int(fields[1]))], contact, n_contacts, n_contacts/(nl+1))
			ppv_score.append(n_contacts/(nl+1))

	ppv_x = [x for x in range(1, len(ppv_score)+1)]
	plt.subplot(111)
	plt.plot(ppv_x, ppv_score)
	plt.show()
