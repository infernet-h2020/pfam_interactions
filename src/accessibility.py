from support import *

def compute_accessibility(pdbname, pdb_path, pdb_dca_resids, allowed_residues, vmd_path, results_folder):
	print("SASA:")
	residues_linear = []
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbname, pdb_path)
	distance_threshold = 20
	for model in structure:
		for chain1 in allowed_residues:
			c1 = model[chain1]
			for r1 in c1:
				resid1 = r1.id[1]
				if resid1 not in allowed_residues[chain1]:
					continue
				residues_linear.append((chain1, resid1))

	N = len(residues_linear)

#	print(residues_linear)
#	print(allowed_residues)

	SASA = {}
	for i1 in range(N):
		c1, r1 = residues_linear[i1]
		for pfam_acc_ann1, dca_i1 in pdb_dca_resids[(c1, r1)]:
#			if dca_i1 in accessibility_dca_resids:
			pfam_acc1, cmult1 = pfam_acc_ann1.split('_')
			os.system("sed 's/XXX/chain {0} and resid {1}/' support_files/fixed_sasa_template.tcl > support_files/fixed_sasa.tcl".format(c1, r1))
			os.system(vmd_path + " " + pdb_path + " -dispdev text -e support_files/fixed_sasa.tcl &> /dev/null")
			with open("SASA_output.dat") as output_file:
				for line in output_file:
					if not dca_i1 in SASA:
						SASA[dca_i1] = []
					SASA[dca_i1].append((pfam_acc_ann1, pdbname+' '+c1+' '+str(r1), float(line.strip())))
					break
	os.remove("SASA_output.dat")

	acessibility_filename = results_folder + "{0}_SASA.txt".format(pdbname)
	with open(acessibility_filename, 'w') as accessibility_file:
		for dca_i in SASA:
			for x in SASA[dca_i]:
				accessibility_file.write("{0}\t{1}\t{2}\t{3}\n".format(dca_i, x[0], x[1], x[2]))
	print(acessibility_filename)
	return SASA
