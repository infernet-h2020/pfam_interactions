from support import *

def mindistance(mind_pdbs, inpfam1, inpfam2, restrict_comparison, results_folder, dca_filename, pdbmap_filename, main_backmap_table, with_offset=True, percdist=False):
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
	totlength = int(text[-1].split()[0]) + max([int(x.split()[1]) for x in text]) - min([int(x.split()[1]) for x in text]) + 1

	recs = []
	if dca_filename and os.path.exists(dca_filename):
		with open(dca_filename) as dca_file:
			for line in dca_file:
				if not line or line.strip().startswith("#"):
					continue
				fields = line.split()
				recs.append((int(fields[0]), int(fields[1]), float(fields[2])))
#				print(dca_filename, int(fields[0]), int(fields[1]), float(fields[2]))
		recs = sorted(recs, key= lambda x: -x[2])
	else:
		for i in range(1, totlength+1):
			for j in range(i+1, totlength+1):
				recs.append((i, j, ""))

	filename_lists = []
	descrs = []
	lcount = -1
	if restrict_comparison:
		groupdef_file = open(results_folder + '/group_definitions.txt', 'w')
		with open(restrict_comparison) as restrict_comparison_file:
			for line in restrict_comparison_file:
				if not line:
					continue

				flist = []
				fields = line.split()
				descrs.append(line.strip())
				lcount += 1
				groupdef_file.write("group {0}: {1}".format(lcount, line))
				if fields[0] == "self":
					if with_offset:
						print("Error (logic: no self in pf1 != pf2 comparisons): in {0}, line\n{1}".format(restrict_comparison, line))
						exit(1)
					if len(fields) > 1:
						textpair = re.findall("\[(.*?)\]", line)
						if len(textpair) != 1:
							print("Error (syntax: only one group allowed): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						group = [x.strip() for x in textpair[0].split(',')]
						for fname in filenames_list:
							if fname.split('_')[-3][0] in group and fname.split('_')[-3] == fname.split('_')[-2]:
								flist.append(fname)
					else:
						for fname in filenames_list:
							if fname.split('_')[-3] == fname.split('_')[-2]:
								flist.append(fname)
				elif fields[0] == "intra" or (fields[0] == "intra-not-self" and with_offset):
					if len(fields) > 1:
						textpair = re.findall("\[(.*?)\]", line)
						if len(textpair) != 1:
							print("Error (syntax: only one group allowed): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						group = [x.strip() for x in textpair[0].split(',')]
						for fname in filenames_list:
							if fname.split('_')[-3][0] in group and fname.split('_')[-3][0] == fname.split('_')[-2][0]:
								flist.append(fname)
					else:
						for fname in filenames_list:
							if fname.split('_')[-3][0] == fname.split('_')[-2][0]:
								flist.append(fname)
				elif fields[0] == "intra-not-self" and not with_offset:
					if len(fields) > 1:
						textpair = re.findall("\[(.*?)\]", line)
						if len(textpair) != 1:
							print("Error (syntax: only one group allowed): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						group = [x.strip() for x in textpair[0].split(',')]
						for fname in filenames_list:
							if fname.split('_')[-3][0] in group and fname.split('_')[-3][0] == fname.split('_')[-2][0] and fname.split('_')[-3][1:] != fname.split('_')[-2][1:]:
								flist.append(fname)
					else:
						for fname in filenames_list:
							if fname.split('_')[-3][0] == fname.split('_')[-2][0] and fname.split('_')[-3][1:] != fname.split('_')[-2][1:]:
								flist.append(fname)
				elif fields[0] == "inter":
					if len(fields) > 1:
						textpair = re.findall("\[(.*?)\]", line)
						if len(textpair) != 1:
							print("Error (syntax: only one group allowed): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						group = [x.strip() for x in textpair[0].split(',')]
						for fname in filenames_list:
							if fname.split('_')[-3][0] in group and fname.split('_')[-3][0] != fname.split('_')[-2][0]:
								flist.append(fname)
					else:
						for fname in filenames_list:
							if fname.split('_')[-3][0] != fname.split('_')[-2][0]:
								flist.append(fname)
				elif fields[0] == "not-self":
					if len(fields) > 1:
						textpair = re.findall("\[(.*?)\]", line)
						if len(textpair) != 1:
							print("Error (syntax: only one group allowed): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						group = [x.strip() for x in textpair[0].split(',')]
						for fname in filenames_list:
							if fname.split('_')[-3][0] in group and fname.split('_')[-3] != fname.split('_')[-2]:
								flist.append(fname)
					else:
						for fname in filenames_list:
							if fname.split('_')[-3] != fname.split('_')[-2]:
								flist.append(fname)
				elif fields[0] == "inter-groups":
					if len(fields) > 1:
						textpairs = re.findall("\[(.*?)\]", line)
						if len(textpairs) != 2:
							print("Error (syntax: two groups required): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						group1 = [x.strip() for x in textpairs[0].split(',')]
						group2 = [x.strip() for x in textpairs[1].split(',')]
						if set(group1) & set(group2):
							print("Error (syntax: intra groups should not intersect): in {0}, line\n{1}".format(restrict_comparison, line))
							exit(1)
						for fname in filenames_list:
							print(fname)
							if (fname.split('_')[-3][0] in group1 and fname.split('_')[-2][0] in group2) or (fname.split('_')[-2][0] in group2 and fname.split('_')[-3][0] in group1):
								flist.append(fname)
					else:
						print("Error (syntax: two groups required): in {0}, line\n{1}".format(restrict_comparison, line))
						exit(1)

				if not flist:
					print("Error (no match found): in {0}, line\n{1}".format(restrict_comparison, line))
					exit(1)

				filename_lists.append(flist)
				flist = []

			if not filename_lists:
				print("Error (empty file): in {0}".format(restrict_comparison))
				exit(1)
			groupdef_file.close()
	else:
		with open(results_folder + '/group_definitions.txt', 'w') as gf:
			gf.write('group 0: all-vs-all\n')
		filename_lists.append(filenames_list)
		descrs.append('all-vs-all')


	for fli, flist in enumerate(filename_lists):
		if percdist:
			if dca_filename and os.path.exists(dca_filename):
				out_filename_3 = results_folder + "/{0}_{1}_g{2}_correlation_percdist_wDCA.txt".format(inpfam1, inpfam2, fli)
				out_filename_4 = results_folder + "/{0}_{1}_g{2}_correlation_percdist_wDCA_fmt.txt".format(inpfam1, inpfam2, fli)
#				print("\n\nStrongest DCA signals mapped on average pair distance, group {0} ({1}):".format(fli, descrs[fli]))
				sorted_dca_filename = os.path.dirname(dca_filename) + '/.sorted_' + os.path.basename(dca_filename)
			else:
				out_filename_3 = results_folder + "/{0}_{1}_g{2}_correlation_percdist.txt".format(inpfam1, inpfam2, fli)
#				print("\n\nFirst DCA paired indices mapped on average pair distance, group {0} ({1}):".format(fli, descrs[fli]))
			out_filename_fig_2 = results_folder + "/{0}_{1}_g{2}_correlation_percdist.png".format(inpfam1, inpfam2, fli)
		if dca_filename and os.path.exists(dca_filename):
			out_filename = results_folder + "/{0}_{1}_g{2}_correlation_mindist_wDCA.txt".format(inpfam1, inpfam2, fli)
			out_filename_2 = results_folder + "/{0}_{1}_g{2}_correlation_mindist_wDCA_fmt.txt".format(inpfam1, inpfam2, fli)
			print("\n\nStrongest DCA signals mapped on closest residue pairs, group {0} ({1}):".format(fli, descrs[fli]))
			sorted_dca_filename = os.path.dirname(dca_filename) + '/.sorted_' + os.path.basename(dca_filename)
		else:
			out_filename = results_folder + "/{0}_{1}_g{2}_correlation_mindist.txt".format(inpfam1, inpfam2, fli)
			print("\n\nFirst DCA paired indices mapped on closest residue pairs, group {0} ({1}):".format(fli, descrs[fli]))
		out_filename_fig = results_folder + "/{0}_{1}_g{2}_correlation_mindist.png".format(inpfam1, inpfam2, fli)

		distcomp = np.ones((len(flist), totlength, totlength))*10000
		rdict = {}
		rdict2 = {}
		for nf, res_filename in enumerate(flist):
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
		if percdist:
			if with_offset:
				distfig = np.ones((offset, distcomp.shape[1]-offset))/2
			else:
				distfig = np.ones(distcomp.shape[1:])/2
			if dca_filename and os.path.exists(dca_filename):
				out_file_4 = open(out_filename_4, 'w')
			with open(out_filename_3, 'w') as out_file:
				for nl, v in enumerate(recs):
					ind1, ind2, score = v
					if with_offset:
						ind1a = ind1 - 1
						ind2a = ind2 - offset - 1
						if not (ind1 <= offset and ind2 > offset):
							continue
					else:
						ind1a = ind1 - 1
						ind2a = ind2 - 1
					if np.min(distcomp[:, ind1-1, ind2-1]) == 10000:
						perc = -1
						distfig[ind1a, ind2a] = 0.25
					else:
						sigmoid_dist = [linear_response(x) for x in distcomp[:, ind1-1, ind2-1]]
						perc = sum(sigmoid_dist)/len(sigmoid_dist)
						distfig[ind1a, ind2a] = 0.5 + perc/2
#					print(with_offset, offset, ind1, ind2, ind1a, ind2a, distfig[ind1a, ind2a])
					if perc >= 0:
						i += 1
						if score == "":
							optarg = ""
						else:
							optarg = "{0: .5f}\t{1}".format(score, nl)
						if i < imax:
							if i == 1:
								out_file.write("# Files can be retrieved in {0}\n".format(os.path.dirname(flist[0])))
						out_file.write("{0}\t{1}\t{2}\t{3}\n".format(ind1, ind2, perc, optarg))
						if dca_filename and os.path.exists(dca_filename):
							out_file_4.write("{0}\t{1}\t{2}\t{3}\n".format(ind1, ind2, score, perc))

#			print(distfig)
#			with open("d.txt", 'w') as f:
#				for i in range(distfig.shape[0]):
#					for j in range(distfig.shape[1]):
#						f.write("{0}\t{1}\t{2}\n".format(i, j, distfig[i,j]))
			plt.imshow(distfig, cmap='RdGy_r', interpolation='none', vmin=0, vmax=1)
			plt.savefig(out_filename_fig_2)

		i = 0	
		if dca_filename and os.path.exists(dca_filename):
			out_file_2 = open(out_filename_2, 'w')
		with open(out_filename, 'w') as out_file:
			if with_offset:
				distfig = np.ones((offset, distcomp.shape[1]-offset))/2
			else:
				distfig = np.ones(distcomp.shape[1:])/2
			for nl, v in enumerate(recs):
				ind1, ind2, score = v
				if with_offset:
					ind1a = ind1 - 1
					ind2a = ind2 - offset - 1
					if not (ind1 <= offset and ind2 > offset):
						continue
				else:
					ind1a = ind1 - 1
					ind2a = ind2 - 1
				if np.min(distcomp[:, ind1-1, ind2-1]) == 10000:
					argmind = ""
					argmind2 = ""
					distfig[ind1a, ind2a] = 0.25
				else:
					nfmin = np.argmin(distcomp[:, ind1-1, ind2-1])
					argmind = (flist[nfmin], rdict[(flist[nfmin], ind1, ind2)], score)
					argmind2 = (ind1, ind2, score, rdict2[(flist[nfmin], ind1, ind2)])
					distfig[ind1a, ind2a] = 0.5 + linear_response(np.min(distcomp[:, ind1-1, ind2-1]))

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
		print("\n\n")

		plt.imshow(distfig, cmap='RdGy_r', interpolation='none', vmin=0, vmax=1)
		plt.savefig(out_filename_fig)
