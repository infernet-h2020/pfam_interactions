# Complete domain-domain interaction database files

The files contained in this folder summarize the Pfam domain-domain interaction
benchmark. They are based on the information found on [PDBfam][pdbfam], which reports
the location of every Pfam domain in each PDB structure.
What we are really interested in is not just the presence of two Pfam domains in a
PDB structure or PDB chain, but whether the regions where the Pfams are mapped are
in contact or not. Our original contribution was thus to perform this analysis on 
the whole PDB.

The table `PfamPfam_contacts_through_all_PDBs_table.txt` summarizes all the important
information on this subject.
It is organized in order to be queried with a pair of Pfam domains, which can also
coincide in case we are intersted in same-domain inter-domain contacts (such as in 
the case of homodimers).

The first two columns report the pair of Pfam domains that are being queried.
The third column reports the name of each PDB structure containing the two Pfam domains,
and the fourth and the fifth contain the chain names where the first and second (respectively)
Pfam domain has been found. 
The sixth column tells whether the two instances of the two domains are in contact. The
keyword "Contact" means there are at least 5 pairs of residues belonging to the two domains
that have a distance < 5 angstrom.
The seventh and eighth columns report the number of instances of the first and second
Pfam domain (respectively) in the chain that is being analyzed. This is important to
understand whether the observed contacts are being biased by other occurrences of the
same domain nearby. 
Columns ninth and tenth tell whether in the chain of the first Pfam domain can be found
at least an instance of the second Pfam domain, and viceversa.
Columns 11 to 13 report, if present, an example of PDB structure where the first Pfam
domain is present in absence of the second Pfam domain, and columns 14 to 16 report
the opposite.
Lastly, the final two columns report the number of sequences contained in the full Pfam
alignment corresponding to the two Pfam domains.

[pdbfam]: http://dunbrack2.fccc.edu/ProtCid/PDBfam/
