from support import *
import freesasa

def compute_accessibility(pdbname, pdb_path, pdb_dca_resids, backmap_table, allowed_residues, vmd_path, results_folder, accessibilities_by_domain):
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
	chains = [x.get_id() for x in structure[0].get_chains()]
	residues = {}
	for c in chains:
		residues[c] = [int(x.get_id()[1]) for x in structure[0][c].get_residues()]
	distance_threshold = 20

	if accessibilities_by_domain:
		for pfam_acc, mapping, unip in backmap_table:
			loc_coords = []
			loc_names = []
			loc_radii = []
			for c, r in [x for x in pdb_dca_resids if x[0]==mapping[0] and x[1]>=mapping[1][0] and x[1]<=mapping[1][1]]:
				if (c not in chains) or (r not in residues[c]):
					continue
				res = structure[0][c][r]
				resname = res.get_resname()
				for atom in res:
					atname = atom.get_name()[0]
					if atname not in radii_d:
						continue
					loc_names.append((c, (r, resname), atname))
					loc_coords.append(atom.get_coord())
					loc_radii.append(radii_d[atname])
			coords.append(loc_coords)
			names.append(loc_names)
			radii.append(loc_radii)
	else:
		loc_coords = []
		loc_names = []
		loc_radii = []
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
						loc_names.append((chain1, (resid1, resname), atname1))
						loc_coords.append(a1.get_coord())
						loc_radii.append(radii_d[atname1])
		coords.append(loc_coords)
		names.append(loc_names)
		radii.append(loc_radii)

	sasas = []
	for i in range(len(coords)):
		loc_coords = coords[i]
		loc_names = names[i]
		loc_radii = radii[i]

		N = len(residues_linear)
		loc_coords = np.array(loc_coords).flatten()
		loc_radii = np.array(loc_radii)

#		print(loc_coords)	
		Results = freesasa.calcCoord(loc_coords, loc_radii)
		sasas.append([Results.atomArea(i) for i in range(Results.nAtoms())])


#	print(names)
	SASA = {}
	if accessibilities_by_domain:
		for ie, entry in enumerate(backmap_table):
			pfam_acc, mapping, unip = entry
			loc_names = names[ie]
			SASA[pfam_acc] = {}
			for i, name in enumerate(loc_names):
				ch, res, at = name
				r, rname = res
				idx = [i for i, x in enumerate(pdb_dca_resids[(ch, r)]) if x[0] == pfam_acc][0]
				dca_i = pdb_dca_resids[(ch, r)][idx][1]
				if (dca_i, ch, r) not in SASA[pfam_acc]:
					SASA[pfam_acc][(dca_i, ch, r)] = [0, 0]
				SASA[pfam_acc][(dca_i, ch, r)][0] += sasas[ie][i]
				rname1 = from3to1(rname)
				if rname1 not in alaxala:
					SASA[pfam_acc][(dca_i, ch, r)][1] = 0
				else:
					SASA[pfam_acc][(dca_i, ch, r)][1] += sasas[ie][i]/alaxala[rname1]
	else:
		loc_names = names[0]
		for ie, entry in enumerate(backmap_table):
			pfam_acc, mapping, unip = entry
			SASA[pfam_acc] = {}
			for i, name in enumerate(loc_names):
				ch, res, at = name
				r, rname = res
				idx_list = [i for i, x in enumerate(pdb_dca_resids[(ch, r)]) if x[0] == pfam_acc]
				if idx_list:
					idx = idx_list[0]
				else:
					continue
				dca_i = pdb_dca_resids[(ch, r)][idx][1]
				if not (ch == mapping[0] and r >= mapping[1][0] and r <= mapping[1][1]):
					continue
				if (dca_i, ch, r) not in SASA[pfam_acc]:
					SASA[pfam_acc][(dca_i, ch, r)] = [0, 0]
				SASA[pfam_acc][(dca_i, ch, r)][0] += sasas[0][i]
				rname1 = from3to1(rname)
				if rname1 not in alaxala:
					SASA[pfam_acc][(dca_i, ch, r)][1] = 0
				else:
#					print(dca_i, ch, r, rname, at, sasas[0][i]/alaxala[rname1])
					SASA[pfam_acc][(dca_i, ch, r)][1] += sasas[0][i]/alaxala[rname1]
				
	acessibility_filename = results_folder + "{0}_SASA.txt".format(pdbname)
	with open(acessibility_filename, 'w') as accessibility_file:
		for ie, entry in enumerate(backmap_table):
			pfam_acc, mapping, unip = entry
			for dca_i, ch, r in sorted(list(SASA[pfam_acc].keys()), key= lambda x : x[0]):
				if SASA[pfam_acc][(dca_i, ch, r)][1] > 1:
					SASA[pfam_acc][(dca_i, ch, r)][1] = 1
				accessibility_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(pfam_acc, dca_i, ch, r, SASA[pfam_acc][(dca_i, ch, r)][0], SASA[pfam_acc][(dca_i, ch, r)][1]))

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
