import os
import sys
import uniprotACC_to_uniprotID

def stockholm_to_fasta(stockholm_file_list, fasta_output_folder):
	if not os.path.exists(fasta_output_folder):
		print("ERROR: fasta output folder {0} does not exist".format(fasta_output_folder))
		exit(1)

	for pf_filename in stockholm_file_list:
		if not os.path.exists(pf_filename):
			print("ERROR: stockholm file {0} does not exist".format(pf_filename))
			exit(1)
		with open(pf_filename) as pf_file:
			query_list = []
			for pfline in pf_file:
				if not pfline or pfline.startswith("#") or pfline.startswith('//'):
					continue
				pffields = pfline.split()
				uniname = pffields[0].split('.')[0]
				query_list.append(uniname)

		with open(pf_filename) as pf_file:
			pfout_filename = fasta_output_folder + "".join([x+"." for x in os.path.basename(pf_filename).split(".")[:-1]])[:-1] + ".fasta"
			with open(pfout_filename, 'w') as pfout_file:
				n = 0
				t = 0
				c = uniprotACC_to_uniprotID.uniprotACC_to_uniprotID(query_list)
				for pfline in pf_file:
					if not pfline or pfline.startswith("#") or pfline.startswith('//'):
						continue
					t += 1
					pffields = pfline.split()
					uniname = pffields[0].split('.')[0]
					seq = pffields[1]
					if uniname in c:
						uniprot_name_complete = c[uniname]
						species = uniprot_name_complete.split("_")[0]
						pfout_file.write(">tr|{0}|{1}\n{2}\n".format(species, uniprot_name_complete, seq))
					else:
						n += 1
				if t != 0:
					print("{0} sequences not found: {1} out of {2} {3:6.3}%".format(pf_filename, n, t, 100*n/t))
				else:
					print("{0} WARNING: 0 sequences found".format(pf_filename))


if __name__ == "__main__":
	stockholm_filename = sys.argv[1]
	stockholm_file_list = [stockholm_filename]
	stockholm_to_fasta(stockholm_file_list, "./")
