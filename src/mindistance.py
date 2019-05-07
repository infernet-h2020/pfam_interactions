from support import *

def mindistance(mind_pdbs, inpfam1, inpfam2, only_intra, only_inter, results_folder, dca_filename, pfam_pfam_filename, with_offset=True):
	text = []
	for pdbname in mind_pdbs:
#		print("grep {0} {1} | grep {2} | grep {3} | awk '{{if ($1==\"{0}\") {{print \"{3}\", $4, $5}} else {{print \"{3}\", $5, $4}} }}'".format(inpfam1, pfam_pfam_filename, inpfam2, pdbname))
		text += subprocess.run(["grep {0} {1} | grep {2} | grep {3} | awk '{{if ($1==\"{0}\") {{print \"{3}\", $4, $5}} else {{print \"{3}\", $5, $4}} }}'".format(inpfam1, pfam_pfam_filename, inpfam2, pdbname)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
#		subprocess.run(["grep {0} {1} | grep {2} | grep {3}".format(inpfam1, pfam_pfam_filename, inpfam2, pdbname)], shell=True)

	pairs = []
	for line in text:
		if not line:
			continue
		pdb, p1, p2 = line.strip().split()
		if only_intra and p1 == p2:
			pairs.append((pdb, p1, p2))
		elif only_inter and p1 != p2:
			pairs.append((pdb, p1, p2))
		elif not (only_intra or only_inter):
			pairs.append((pdb, p1, p2))

	filenames_list = []
	for pair in pairs:
		fname = results_folder + "/{0}_{1}_{2}_{3}*_{4}*.txt".format(inpfam1, inpfam2, pair[0].lower(), pair[1], pair[2])
		text = subprocess.run(["ls", fname], stderr=open("/dev/null", 'w'), stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
		for line in text:
			if not line:
				continue
			filenames_list.append(line.strip())
	
#	print("FILES", filenames_list)

	if with_offset:
		with open(filenames_list[0]) as one_file:
			for line in one_file:
				if not line:
					continue
				offset = int(line.split()[0])

	imax = 50
	i = 0
	out_filename = results_folder + "/{0}_{1}_correlation.txt".format(inpfam1, inpfam2)
	if only_intra:
		warner = " (only intra-chain)"
	elif only_inter:
		warner = " (only inter-chain)"
	else:
		warner = ""
	print("\n\nStrongest DCA signals mapped on closest residue pairs{0}:".format(warner))
	sorted_dca_filename = os.path.dirname(dca_filename) + '/.sorted_' + os.path.basename(dca_filename)
	with open(dca_filename) as dca_file:
		recs = []
		for line in dca_file:
			if not line or line.strip().startswith("#"):
				continue
			fields = line.split()
			recs.append((int(fields[0]), int(fields[1]), float(fields[2])))

	with open(out_filename, 'w') as out_file:
		for nl, v in enumerate(sorted(recs, key= lambda x: -x[2])):
			ind1, ind2, score = v
			if with_offset and (not (ind1 <= offset and ind2 > offset)):
				continue
			outs = []
			for res_filename in filenames_list:
				text = subprocess.run(["grep '^\s*{0}\s\s*{1}\s' {2}".format(ind1, ind2, res_filename)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
				if text:
					text = text[:-1]
					outs.append((res_filename, text, score))
			mind = 1000000
			argmind = ""
			for o in outs:
				d = o[1].split()[6]
				if d == "None":
					continue
				elif float(d) < mind:
					mind = float(d)
					argmind = o
			if argmind:
				i += 1
				if i < imax:
					print("{0}\t{1}\t{2: .5f}\t{3}".format(os.path.basename(argmind[0]), argmind[1], argmind[2], nl))
					if i == 1:
						out_file.write("# Files can be retrieved in {0}\n".format(os.path.dirname(argmind[0])))
				out_file.write("{0}\t{1}\t{2: .5f}\t{3}\n".format(os.path.basename(argmind[0]), argmind[1], argmind[2], nl))

	print("\nFind this table in {0}".format(out_filename))
