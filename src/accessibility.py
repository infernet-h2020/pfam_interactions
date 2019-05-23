from support import *
import freesasa

def compute_accessibility(pdbname, pdb_path, pdb_dca_resids, allowed_residues, vmd_path, results_folder):
	print("SASA:")

	radii_d = {}
	this_path = os.path.dirname(os.path.abspath(__file__)) + '/'
	with open(this_path + 'support_files/radii.txt') as radii_file:
		for line in radii_file:
			fields = line.split()
			radii_d[fields[0]] = float(fields[1])

	alaxala = {}
	with open(this_path + 'support_files/alaxala_sasa.txt') as alaxala_file:
		for line in alaxala_file:
			fields = line.split()
			alaxala[fields[0]] = float(fields[1])

	residues_linear = []
	coords = []
	names = []
	radii = []
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbname, pdb_path)
	distance_threshold = 20
	for model in structure:
		for chain1 in allowed_residues:
			c1 = model[chain1]
			for r1 in c1:
				resid1 = r1.id[1]
				resname = r1.get_resname()
				if resid1 not in allowed_residues[chain1]:
					continue
				if (chain1, resid1) not in residues_linear:
					residues_linear.append((chain1, resid1))
				for a1 in r1:
					atname1 = a1.get_name()[0]
					if atname1 not in radii_d:
						continue
					names.append((chain1, (resid1, resname), atname1))
					coords.append(a1.get_coord())
					radii.append(radii_d[atname1])
	N = len(residues_linear)
	coords = np.array(coords).flatten()
	radii = np.array(radii)

	Results = freesasa.calcCoord(coords, radii)
	sasas = [Results.atomArea(i) for i in range(Results.nAtoms())]
	pre_SASA = {}
	pre_rSASA = {}
	for i, name  in enumerate(names):
		ch, res, at = name
		r, rname = res
		if (ch, r) not in pre_SASA:
			pre_SASA[(ch, r)] = 0
			pre_rSASA[(ch, r)] = 0
		pre_SASA[(ch, r)] += sasas[i]
		rname1 = from3to1(rname)
		if rname1 not in alaxala:
			pre_rSASA[(ch, r)] = 0
		else:
			pre_rSASA[(ch, r)] += sasas[i]/alaxala[rname1]


	SASA = {}
	for ch, r in pre_SASA:
		for pfam_acc_ann1, dca_i1 in pdb_dca_resids[(ch, r)]:
			pfam_acc1, cmult1 = pfam_acc_ann1.split('_')
			if not dca_i1 in SASA:
				SASA[dca_i1] = []
			SASA[dca_i1].append((pfam_acc_ann1, pdbname+' '+ch+' '+str(r), pre_SASA[(ch, r)], min(1, pre_rSASA[(ch, r)])))

	accessibility = []
	acessibility_filename = results_folder + "{0}_SASA.txt".format(pdbname)
	with open(acessibility_filename, 'w') as accessibility_file:
		for dca_i in SASA:
			for x in SASA[dca_i]:
				accessibility.append(("{0}\t{1}\t{2}\t{3}\t{4}\n".format(dca_i, x[0], x[1], x[2], x[3]), x, dca_i))
		for entry in sorted(accessibility, key= lambda x: (x[1], x[2])):
			accessibility_file.write(entry[0])
	print(acessibility_filename)

	return SASA

#	print(residues_linear)
#	print(allowed_residues)

	if not os.path.exists(vmd_path):
		print("CRITICAL: path {0} not found.\n\tIf it is an alias, please expand it.".format(vmd_path))
		exit(1)
	else:
		vmd_path = vmd_path.replace(" ", "\ ")

	SASA = {}
	for i1 in range(N):
		c1, r1 = residues_linear[i1]
		print(pdb_dca_resids[(c1, r1)])
		for pfam_acc_ann1, dca_i1 in pdb_dca_resids[(c1, r1)]:
#			if dca_i1 in accessibility_dca_resids:
			pfam_acc1, cmult1 = pfam_acc_ann1.split('_')
			os.system("sed 's/XXX/chain {0} and resid {1}/' {2}/support_files/fixed_sasa_template.tcl | sed 's=YYY={3}=g' > {2}/support_files/fixed_sasa.tcl".format(c1, r1, this_path, results_folder))
			sp = subprocess.run([vmd_path + " " + pdb_path + " -dispdev text -e {0}/support_files/fixed_sasa.tcl > /dev/null 2>&1".format(this_path)], shell=True, executable='/bin/bash')
			with open(results_folder + "SASA_output.dat") as output_file:
				for line in output_file:
					if not dca_i1 in SASA:
						SASA[dca_i1] = []
					SASA[dca_i1].append((pfam_acc_ann1, pdbname+' '+c1+' '+str(r1), float(line.strip())))
					break
	os.remove(results_folder + "SASA_output.dat")

	acessibility_filename = results_folder + "{0}_SASA.txt".format(pdbname)
	with open(acessibility_filename, 'w') as accessibility_file:
		for dca_i in SASA:
			for x in SASA[dca_i]:
				accessibility_file.write("{0}\t{1}\t{2}\t{3}\n".format(dca_i, x[0], x[1], x[2]))
	print(acessibility_filename)
	return SASA
