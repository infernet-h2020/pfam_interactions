from support import *

def mindistance(mind_pdbs, inpfam1, inpfam2, only_intra, only_inter, results_folder, dca_filename, pdbmap_filename, main_backmap_table, with_offset=True):
	pairs = []
	for pdbname in mind_pdbs:
		chains = {}
		for inpfam in [inpfam1, inpfam2]:
			chains[inpfam] = subprocess.run(["zgrep {0} {1} | grep {2} | awk '{{print substr($2,1,1)}}'".format(pdbname.upper(), pdbmap_filename, inpfam)], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
			chains[inpfam] = [x for x in chains[inpfam] if x.strip()]
		for ch1 in chains[inpfam1]:
			for ch2 in chains[inpfam2]:
				pairs.append((pdbname, ch1, ch2))

	filenames_list = []
	for pair in pairs:
		fname = results_folder + "/{0}_{1}_{2}_{3}*_{4}*_*.txt".format(inpfam1, inpfam2, pair[0].lower(), pair[1], pair[2])
		text = subprocess.run(["ls {0}".format(fname)], stderr=open("/dev/null", 'w'), stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')
		for line in text:
			if not line:
				continue
			filenames_list.append(line.strip())
	
	one_file = open(filenames_list[0])
	text = [l for l in one_file.readlines() if l.strip()]
	if with_offset:
		offset = int(text[-1].split()[0])
	totlength = int(text[-1].split()[1])

	if only_intra:
		warner = " (only intra-chain)"
	elif only_inter:
		warner = " (only inter-chain)"
	else:
		warner = ""

	recs = []
	if dca_filename and os.path.exists(dca_filename):
		out_filename = results_folder + "/{0}_{1}_correlation_wDCA.txt".format(inpfam1, inpfam2)
		out_filename_2 = results_folder + "/{0}_{1}_correlation_wDCA_fmt.txt".format(inpfam1, inpfam2)
		print("\n\nStrongest DCA signals mapped on closest residue pairs{0}:".format(warner))
		sorted_dca_filename = os.path.dirname(dca_filename) + '/.sorted_' + os.path.basename(dca_filename)
		with open(dca_filename) as dca_file:
			for line in dca_file:
				if not line or line.strip().startswith("#"):
					continue
				fields = line.split()
				recs.append((int(fields[0]), int(fields[1]), float(fields[2])))
#				print(dca_filename, int(fields[0]), int(fields[1]), float(fields[2]))
		recs = sorted(recs, key= lambda x: -x[2])
	else:
		out_filename = results_folder + "/{0}_{1}_correlation.txt".format(inpfam1, inpfam2)
		print("\n\nFirst DCA paired indices mapped on closest residue pairs{0}:".format(warner))
		for i in range(1, totlength+1):
			for j in range(i+1, totlength+1):
				recs.append((i, j, ""))

	distcomp = np.ones((len(filenames_list), totlength, totlength))*10000
	rdict = {}
	rdict2 = {}
	for nf, res_filename in enumerate(filenames_list):
		with open(res_filename) as res_file:
			for line in res_file:
				fields = line.split()
				if fields[6] != "None":
					if fields[6] == "MAX":
						distcomp[nf, int(fields[0])-1, int(fields[1])-1] = float(9999)
					else:
						distcomp[nf, int(fields[0])-1, int(fields[1])-1] = float(fields[6])
					rdict[(res_filename, int(fields[0]), int(fields[1]))] = line.strip()
					pfam_acc1, pfam_acc2, pdbname, chain_ann1, chain_ann2, _ = os.path.basename(res_filename).split("_")
					pfam_acc_ann1, pfam_acc_ann2 = pfam_acc1 + "_" + chain_ann1, pfam_acc2 + "_" + chain_ann2
					init1, end1 = [x[1][1] for x in main_backmap_table[pdbname] if x[0] == pfam_acc_ann1][0]
					init2, end2 = [x[1][1] for x in main_backmap_table[pdbname] if x[0] == pfam_acc_ann2][0]
					if fields[6] == "MAX":
						rdict2[(res_filename, int(fields[0]), int(fields[1]))] = ("MAX", pdbname+"_"+chain_ann1[0]+":"+str(init1)+"-"+str(end1)+"_"+chain_ann2[0]+":"+str(init2)+"-"+str(end2), fields[3], fields[5])
					else:
						rdict2[(res_filename, int(fields[0]), int(fields[1]))] = (fields[6], pdbname+"_"+chain_ann1[0]+":"+str(init1)+"-"+str(end1)+"_"+chain_ann2[0]+":"+str(init2)+"-"+str(end2), fields[3], fields[5])


	i = 0
	imax = 50
	if dca_filename and os.path.exists(dca_filename):
		out_file_2 = open(out_filename_2, 'w')
	with open(out_filename, 'w') as out_file:
		for nl, v in enumerate(recs):
			ind1, ind2, score = v
			if with_offset and (not (ind1 <= offset and ind2 > offset)):
				continue
			if np.min(distcomp[:, ind1-1, ind2-1]) == 10000:
				argmind = ""
				argmind2 = ""
			else:
				nfmin = np.argmin(distcomp[:, ind1-1, ind2-1])
				argmind = (filenames_list[nfmin], rdict[(filenames_list[nfmin], ind1, ind2)], score)
				argmind2 = (ind1, ind2, score, rdict2[(filenames_list[nfmin], ind1, ind2)])

			if argmind:
				i += 1
				if argmind[2] == "":
					optarg = ""
				else:
					optarg = "{0: .5f}\t{1}".format(argmind[2], nl)
				if i < imax:
					print("{0}\t{1}\t{2}".format(os.path.basename(argmind[0]), argmind[1], optarg))
					if i == 1:
						out_file.write("# Files can be retrieved in {0}\n".format(os.path.dirname(argmind[0])))
				out_file.write("{0}\t{1}\t{2}\n".format(os.path.basename(argmind[0]), argmind[1], optarg))
				if dca_filename and os.path.exists(dca_filename):
					out_file_2.write("{0}\t{1}\t{2}\t{3}\t{4:30}\t{5}\t{6}\n".format(argmind2[0], argmind2[1], argmind2[2], argmind2[3][0], argmind2[3][1], argmind2[3][2], argmind2[3][3]))

	print("\nFind this table in {0}".format(out_filename))
	if dca_filename and os.path.exists(dca_filename):
		out_file_2.close()
		print("\nand {0}".format(out_filename_2))
