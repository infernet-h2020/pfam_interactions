# Pfam domain-domain interaction benchmark

This repository aims to give a complete list of all domain-domain interactions among [Pfam][pfam]
domains present in each [PDB][pdb] structure. The benchmark is complemented by a library of Python 3.6
tools for mapping DCA predictions onto PDB structures.


## Contents 

* complete\_database/ contains the *benchmark tables* mapping each Pfam domain (or pair of Pfam domains)
to the PDB structures containing them

* src/ contains the *dca2pdb* tool suite, which allows a direct mapping of any DCA prediction based on 
Pfam multiple sequence alignments to the chosen PDB structure

* examples/ contains several *examples* for the correct use of all dca2pdb options


[pfam]: https://pfam.xfam.org/
[pdb]: https://www.rcsb.org/
