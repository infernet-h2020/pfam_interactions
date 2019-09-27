from support import *

def calculate_matches(pdbname, inpfam, inpfam1, inpfam2, pdb_pfam_filename):
	print("Matches:")
	pfam_in_pdb = ['', '']
	text = subprocess.run(['zgrep', '{0}'.format(pdbname.upper()), pdb_pfam_filename], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
#	print("zgrep {0} {1}".format(pdbname.upper(), pdb_pfam_filename))
	for line in text:
		if not line:
			continue
		fields = line.split()
		pdbc = fields[2] + '_' + fields[5]
		pfam_acc = fields[3]
		uniprot_acc = fields[4]
		if (inpfam and pfam_acc != inpfam) or (inpfam1 and (pfam_acc != inpfam1 and pfam_acc != inpfam2)) or (pfam_in_pdb[0] and pfam_in_pdb[1]):
			continue
		elif inpfam1 and not pfam_in_pdb[0] and pfam_acc == inpfam1:
			pfam_in_pdb[0] = (pfam_acc, uniprot_acc)
		elif inpfam1 and pfam_acc == inpfam2:
			pfam_in_pdb[1] = (pfam_acc, uniprot_acc)
		elif inpfam:
			pfam_in_pdb[0] = (pfam_acc, uniprot_acc)
		print(line)
#		print(pfam_in_pdb)
	pfam_in_pdb_correct = []
	for p in pfam_in_pdb:
		if p:
			pfam_in_pdb_correct.append(p)
	pfam_in_pdb = pfam_in_pdb_correct
	print("")

	# Consistency checks	
	if (inpfam1 and (inpfam1 not in [x[0] for x in pfam_in_pdb] or inpfam2 not in [x[0] for x in pfam_in_pdb])) or inpfam and (inpfam not in [x[0] for x in pfam_in_pdb]):
		print("ERROR: One or more query Pfams was not found in target PDB")
		exit(1)
	if not pfam_in_pdb:
		print("No Pfams selected!")
		exit(1)

	print("PFAM_IN_PDB", pfam_in_pdb)
	return pfam_in_pdb
