# Pfam domain-domain interaction benchmark

This repository aims to give a complete list of all domain-domain interactions among [Pfam][pfam]
domains present in each [PDB][pdb] structure. The benchmark is complemented by a library of Python 3.6
tools for mapping DCA predictions onto PDB structures.


## Contents 

* db/ contains the **benchmark tables** (interaction\_database/) mapping each Pfam domain or pair of Pfam domains
to the respective PDB structures. It is also the location where helper files from Pfam are downloaded (external\_resources/) 

* src/ contains the **pfam2pdb** tool suite, which allows a direct mapping of any DCA prediction based on 
Pfam's multiple sequence alignments to the chosen PDB structure

* examples/ contains several **examples** for the correct use of all pfam2pdb options (please refer to src/README.md)

* .cache/ is the location where pfam2pdb will store precomputed results to speed up further queries


## Installation

### Requirements

* GCC compiler

The package will decompress the Pfam interaction benchmark tables contained in db/interaction\_database/ and download 
a few large database file from Pfam, storing them in db/external\_resources/. The version of the Pfam database used
is by default **Pfam 32.0**, but this can be easily changed by changing the value of the "version" variable in the
"install.sh" script (note: it must be an integer, versions such as 5.1 are not supported right now).
The version of Pfam can also be changed in a second time with the script "src/change\_pfam\_version.py".

To install, just navigate to the main folder of the package and type

`
./install.sh
`

After the installation, we recommend to add the pfam\_interations/bin/ path to your PATH environment variable by adding 
this line to your ~/.bashrc file (in case you are running bash)

`
export PATH=/your/path/to/pfam\_interations/bin/:$PATH
`

Note: the installation takes much time (2-3 hours), since it has to download ~15 Gb of information from Pfam and then indicize
the files in order to be faster during runtime.


<a name="infernet_logo"/>
<div align="center">
<a href="http://www.infernet.eu/" target="_blank">
<img src="http://www.infernet.eu/wp-content/uploads/2017/03/INFERNET_Wordmark_HR.png" alt="infernet logo" width="200" height="50"></img>
</a>
</div>


This work is supported by [INFERNET](http://www.infernet.eu) : "New algorithms for inference and optimization from large-scale biological data".

<a name="eu_banner"/>
<div align="center">
<a href="https://europa.eu/european-union/index_en" target="_blank">
<img src="http://www.infernet.eu/wp-content/uploads/2017/03/flag_yellow_high.jpg" alt="eu banner" width="40" height="30"></img>
</a>
</div>

<p align="center"><sup>
The INFERNET project is co-funded by the European Union’s H2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement number 734439.
</sup>
</p>






[pfam]: https://pfam.xfam.org/
[pdb]: https://www.rcsb.org/
