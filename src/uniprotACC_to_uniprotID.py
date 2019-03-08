import os
import sys
import urllib
import urllib.parse
import urllib.request
import requests

def uniprotACC_to_uniprotID(query_list):
	url = 'https://www.uniprot.org/uploadlists/'

	result_d = {}

	binning = 10000
	for i in range(int(len(query_list)/binning)+1):
		support_filename = ".uniprot_query.txt"
		support_out_filename = ".uniprot_response.txt"
		with open(support_filename, "w") as support_file:
			for q in query_list[i*binning:(i+1)*binning]:
				support_file.write(q+"\n")

		os.system("perl q.pl {0} > {1}".format(support_filename, support_out_filename))

		with open(support_out_filename) as support_out_file:
			for line in support_out_file:
				if not line:
					continue
				fields = line.split()
				if "_" not in fields[1]:
					continue
				result_d[fields[0]] = fields[1]

		print("Batch done!")
#	print(query_list)
		print(result_d)

	return result_d

if __name__ == "__main__":
	query_filename = sys.argv[1]

	query_list = []
	with open(query_filename, "r") as query_file:
		for line in query_file:
			if not line:
				continue
			query_list.append(line.strip())

	result_d = uniprotACC_to_uniprotID(query_list)

	for k in sorted(list(result_d.keys())):
		print(k, result_d[k])
