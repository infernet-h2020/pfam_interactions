# Pfam domain-domain interaction benchmark

This repository aims to give a complete list of all domain-domain interactions among [Pfam][pfam]
domains present in each [PDB][pdb] structure. The benchmark is complemented by a library of Python 3.6
tools for mapping DCA predictions onto PDB structures.


## Contents 

* db/ contains the **benchmark tables** (interaction\_database/) mapping each Pfam domain or pair of Pfam domains
to the respective PDB structures. 

* src/ contains the **dca2pdb** tool suite, which allows a direct mapping of any DCA prediction based on 
Pfam's multiple sequence alignments to the chosen PDB structure

* examples/ contains several **examples** for the correct use of all dca2pdb options (please refer to src/README.md)


## Installation

The package will decompress the Pfam interaction benchmark tables contained in db/interaction\_database/ and download 
a few large database file from Pfam, storing them in db/external\_resources/. The version of the Pfam database used
is by default **Pfam 32.0**, but this can be easily changed by changing the value of the "version" variable in the
"install.sh" script (note: it must be an integer, versions such as 5.1 are not supported right now).
The version of Pfam can also be changed in a second time with the script "src/change\_pfam\_version.py".

To install, just navigate to the main folder of the package and type

`
./install.sh
`

There are no options to set.

[pfam]: https://pfam.xfam.org/
[pdb]: https://www.rcsb.org/
