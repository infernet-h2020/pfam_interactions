# DCA-2-PDB tools

This suite of tools is aimed to map any DCA prediction based on Pfam multiple sequence 
alignments on the corresponding PDB structure. It allows to:

* check whether two Pfam domains interact in a specific PDB

* test the quality of both single-domain and domain-domain DCA contact predictions

* check whether the residues DCA predicted to make an inter-protein contact are exposed
on the surface of the monomers (especially helpful in case the structure of the
complex is missing)


## Required software

The free molecular visualization software [VMD][vmd] Version 1.9.3 is needed to 
compute the exposed surface of the interacting residues.

## Tutorial - Testing the quality of a DCA contact prediction based on Pfam MSAs

The main scope of the Pfam interaction database presented in this repository is to serve
as a training and test set for any DCA contact prediction algorithm. Hereby we describe
how DCA-2-PDB tools can help tracking and understanding the performance of any DCA
contact prediction algorithm. Although the examples were generated with plmDCA, which was
specifically designed for contact prediction applications, the suite can be used with any
DCA version.
This section describes the main tests that can be done through DCA-2-PDB tools. The suite
does not include the [plmDCA code][plmDCAcode] nor the generation of DCA contact prediction
files, which is in fact one of the expected inputs of the program.

### DCA on a single domain

DCA is often used in order to find contacts within a single protein domain.
As an example, let's take Pfam PF00940 and PDB 3e2e. With the command

`
python3 dca2pdb.py -pf PF02874 -pdb 2wss
`

the backmap algorithm finds all occurrences of the given Pfam family in the given PDB structure,
and for each instance it outputs a text file where all pairs of DCA indices are matched with
their distance in the PDB structure of choice.

The stdout is divided into three parts:

* *Matches* shows all the occurrences of the selected Pfam domain in the selected PDB structure.
In the example, it finds 12 occurrences, one in each of the 12 chains of the structure.

```
Matches:
2146595	321159490	2WSS	PF02874	P00829	E	13	NULL	79	NULL	63	129	2dcf00
2146624	321159490	2WSS	PF02874	P00829	N	13	NULL	79	NULL	63	129	2dcf00
2146636	321159490	2WSS	PF02874	P00829	F	13	NULL	79	NULL	63	129	2dcf00
2146647	321159490	2WSS	PF02874	P00829	D	13	NULL	79	NULL	63	129	2dcf00
2146665	321159490	2WSS	PF02874	P00829	M	13	NULL	79	NULL	63	129	2dcf00
2146708	321159490	2WSS	PF02874	P00829	O	13	NULL	79	NULL	63	129	2dcf00
2146786	321365921	2WSS	PF02874	P19483	B	24	NULL	92	NULL	67	135	2dcf00
2146790	321365921	2WSS	PF02874	P19483	L	30	NULL	92	NULL	73	135	2dcf00
2146832	321365921	2WSS	PF02874	P19483	K	24	NULL	92	NULL	67	135	2dcf00
2146839	321365921	2WSS	PF02874	P19483	A	24	NULL	92	NULL	67	135	2dcf00
2146880	321365921	2WSS	PF02874	P19483	C	27	NULL	92	NULL	70	135	2dcf00
2146885	321365921	2WSS	PF02874	P19483	J	24	NULL	92	NULL	67	135	2dcf00
```

* *Backmap* shows, for each occurrence of the Pfam domain, the PDB residue index range where it
was backmapped. To give the Pfam domain an unambiguous label, the Pfam accession number is
complemented with the 1-letter PDB chain name followed by a progressive index (starting form 1),
should more than one occurrence of the same Pfam domain be found in the same chain.

```
Backmap: 
PF02874_A1: A 26-92
PF02874_B1: B 26-92
PF02874_C1: C 27-92
PF02874_D1: D 13-79
PF02874_E1: E 13-79
PF02874_F1: F 13-79
PF02874_J1: J 26-92
PF02874_K1: K 26-92
PF02874_L1: L 30-92
PF02874_M1: M 13-79
PF02874_N1: N 13-79
PF02874_O1: O 13-79
```

* *Interactions* reports the location of the output interaction file relative to each occurrence
of the selected Pfam domain. The first two fields of the summary table specify which two Pfam 
domains are interacting (this notation will later become more useful), the third field is a
progressive index for the interaction files, the fourth specifies the path where to find the file
and the fifth reports the type of contact prediction to which the file refers to. With only one
Pfam domain, the only type that can be found is "same", which states that the files will relate
the DCA indices of a contact prediction made on a simple, single-domain Pfam MSA to the distances
between residues contained in the same occurrence of the Pfam domain on the selected PDB.

```
Interactions:
	calculating distances	(pickled!)	00:00:00
	preparing support data for calculating interactions	00:00:03
	calculating interactions	00:00:00
PF02874 PF02874 0 PF02874_PF02874_2wss_A1_A1.txt same
PF02874 PF02874 1 PF02874_PF02874_2wss_B1_B1.txt same
PF02874 PF02874 2 PF02874_PF02874_2wss_C1_C1.txt same
PF02874 PF02874 3 PF02874_PF02874_2wss_D1_D1.txt same
PF02874 PF02874 4 PF02874_PF02874_2wss_E1_E1.txt same
PF02874 PF02874 5 PF02874_PF02874_2wss_F1_F1.txt same
PF02874 PF02874 6 PF02874_PF02874_2wss_J1_J1.txt same
PF02874 PF02874 7 PF02874_PF02874_2wss_K1_K1.txt same
PF02874 PF02874 8 PF02874_PF02874_2wss_L1_L1.txt same
PF02874 PF02874 9 PF02874_PF02874_2wss_M1_M1.txt same
PF02874 PF02874 10 PF02874_PF02874_2wss_N1_N1.txt same
PF02874 PF02874 11 PF02874_PF02874_2wss_O1_O1.txt same
```

Each output interaction file reports the following information: DCA index of the two residues considered,
their PDB chain name and residue index, distances (respectively: minimum distance, minimum distance 
between side chains, distance between C\_alphas), their uniprot reference (uniprot accession number, 
residue index and residue type, all condensed in one field and separated by underscores).
This is an extract from "PF02874\_PF02874\_2wss\_A1\_A1.txt":

```
     3       3  A         26    A         26           0.0             0.0             0.0      P19483_69_E     P19483_69_E
     3       4  A         26    A         27         1.322           3.779           3.779      P19483_69_E     P19483_70_E
     3       5  A         26    A         28         3.305           5.582           5.582      P19483_69_E     P19483_71_T
     3       6  A         26    A         29         6.224           7.475           7.475      P19483_69_E     P19483_72_G
     3       7  A         26    A         30         8.901           11.19           11.19      P19483_69_E     P19483_73_R
     3       8  A         26    A         31         12.41           14.19           14.19      P19483_69_E     P19483_74_V
     3       9  A         26    A         32         13.29           14.07           15.66      P19483_69_E     P19483_75_L
     3      10  A         26    A         33         16.62           18.67           18.67      P19483_69_E     P19483_76_S
     3      11  A         26    A         34         20.36           20.36           20.36      P19483_69_E     P19483_77_I
     3      12  A         26    A         35          23.7            23.7            23.7      P19483_69_E     P19483_78_G
     3      13  A         26    A         36         25.69           25.69           25.69      P19483_69_E     P19483_79_D
```

It should be noted that, so far, we haven't made use of any file containing DCA contact predictions.
Indeed, **as long as dca2pdb runs the same version of Pfam used in the DCA predictions**, the program does
not need DCA output files to run.

### True/false contacts and PPV calculation

If a correctly formatted DCA output file is available, we can easily calculate the PPV and other
quality assessment measures. Let's take the example file "examples/plmDCA\_outputs/PF02874\_plmdca.txt"
and the newly calculated "PF02874\_PF02874\_2wss\_A1\_A1.txt":

`
python3 quality_assessment.py -int ../examples/results/PF02874\_PF02874\_2wss\_A1\_A1.txt -dca ../examples/plmDCA\_outputs/PF02874\_plmdca.txt
`

The

### DCA for domain-domain interaction

DCA contact prediction has also been used to infer contacts between pairs of domains. The
multiple sequence alignments relative to the two domains are matched sequence-by-sequence with
the help of a paralog matching algorithm. To generate the example files, we used [ParalogMatching.jl][ParalogMatchingcode].
In `examples/Pfam\_32\_uniprot/plmDCA\_result/`, the file `PF00017_PF00041_plmdca.txt` is a
DCA-formatted output for contact prediction. The prediction was run on the _uniprot_ MSAs of
Pfam Version 32 of the Pfam families PF00017 and PF00041. In a correctly formatted DCA output 
file, columns 1 and 2 contain the DCA index of the column of the Pfam MSA, and column 3 contains 
the DCA score relative to the said residue pair.

In order to test the accuracy of the prediction, we follow the same procedure of the first section
of this tutorial: first, we want to associate the DCA score relative to a pair of residues belonging 
respectively to the first and second Pfam family to the corresponding inter-domain residue-residue 
distance taken from a PDB structure of choice, e.g. 5dc0.
This can be accomplished with

`
python3 dca2pdb.py -pf1 PF00017 -pf2 PF00041 -pdb 5dc0
`

As in the case of single-domain interaction, the program will report on stdout information about
the matches found in the benchmark and the regions of the PDB structures on which the two Pfam 
domains were backmapped. Additionally, it will also show the interactions found, in a table
reporting the Pfam accession numbers of the two domains, a progressive index, the path to the
output file containing the interactions residue-by-residue and the nature of the interaction
("intra" if the two domains are contained in the same chain of the PDB structre, "inter" otherwise).

`
PF00017 PF00041 0 ../examples/results/PF00017_PF00041_5dc0_B1_A1.txt inter
`

The output interaction file reports the following information: DCA index of the two residues considered,
their PDB chain name and residue index, distances (minimum distance, minimum distance between side chains,
distance between C\_alphas), their uniprot reference (uniprot accession number, residue index and residue 
type all condensed in one field).

From the joint MSA obtained by matching the paralog sequences of the two Pfam MSAs of choice, the DCA 
algorithm predicts both contacts internal to each of the two domains and contacts between the two domains.
Yet, our output file only reports those residue pairs where the first residue belongs to the first Pfam domain 
and the second corresponds to the second Pfam domain. For this reason, the first line of the example 
output file PF00017\_PF00041\_5dc0\_B1\_A1.txt reports the information for the residue pair (1, 79), 
the DCA length of Pfam PF00017 being 78.

For calculating the PPV, we can run the usual command:

`
python3 quality_assessment.py -int ../examples/results/PF00017_PF00041_5dc0_B1_A1.txt -dca ../examples/plmDCA_outputs/PF00017_PF00041_plmdca.txt
`

### Combining the coordinates of multiple structures

In many cases, more than one PDB structure containing the chosen pair of Pfam domains is available.
The user can choose to consider any number of PDB structures with the option `-mind`. The option takes
as argument either the path to a text file where a list of PDB codes is specified (one PDB code per
line of the file), or the keyword `all`, which activates a search in the interaction benchmark for
all structures containing the two specified Pfam domains.

As an example, we type:

`
python3 dca2pdb.py -pf1 PF00017 -pf2 PF00041 -dca ../examples/plmDCA_outputs/PF00017_PF00041_plmdca.txt -mind all
`

Note the absence of the `-pdb` keyword, which now is overridden by `-mind`, and the presence of -dca,
which so far was seen only in the quality\_assessment.py algorithm.
The program will tell which PDB structures are considered likely to have both Pfam families:

```
Considering the following PDBs:
3k2m, 3t04, 3uyo, 4je4, 4jeg, 4jmg, 4jmh, 5dc0, 5dc4, 5dc9, 5mtj, 5mtm, 5mtn
```


The program will go through several cycles of backmapping and interaction calculations, one for each 
structure considered. The backmapping routine will fail on some of these structures, and will give the
following error:

`
ERROR: One or more query Pfams was not found in target PDB
`

This is due to the fact that some of the PDB structures that are deemed to contain the two Pfam domains
according to the Pfam interaction benchmark might not in fact strictly contain one or both the domains
according to the internal Pfam database annotations, which adopt stricter definitions. In order to know
more about this topic, please read through the README file of the complete\_databases/ folder. 
At the end of this procedure, the program will inform the user on the structures that have been actually
found to contain the Pfam families:

```
PDBs effectively interpolated:
4je4, 4jeg, 5dc0, 5dc4, 5dc9
```

In case none is found, the algorithm will stop at this point.

When all structures have been searched for interactions, the information of the ones correctly backmapped
is merged, in order to retain for each DCA residue-residue pair the information relative to the PDB structure
where the distance between the two residues is minimal. The justification for this procedure lies in
the signal collected by DCA: if ever two residues make a structurally or functionally important contact,
this could produce coevolutionary signal that DCA could spot. The two residues do not need to be in contact
in any member of the protein family or in any conformation of the same protein, thus it is more correct
to consider the minimum observed distance between any residue-residue pair, should this information be
available.

The algorithm reports in stdout a table of the 50 strongest DCA inter-domain pairs of residues, each one
related to the minimum distance observed in a PDB structure.
Specifically, the table reports: the path to the interaction output filename from where the contact
information was taken, the DCA index of the two residues, their PDB chain and residue index, their 
distances, their uniprot notation, the DCA score of the pair and the rank of the DCA contact.
The last line reports the path where to find this table saved as a text file. 
Here are reported the first entries:

```
Strongest DCA signals mapped on closest residue pairs:
PF00017_PF00041_5dc0_B1_A1.txt	    36	    92	B   	 180	A   	  16	     28.44	     28.44	     28.44	None	P02751_1462_T	0.16359	317
PF00017_PF00041_5dc9_A1_B1.txt	    22	    80	A   	 166	B   	   6	     12.61	     14.96	     15.51	P00519_147_G	P02751_1450_V	0.12042	573
PF00017_PF00041_5dc4_A1_B1.txt	    22	   146	A   	 166	B   	  69	     25.04	     25.04	     25.04	P00519_147_G	P02751_1513_D	0.11715	598
PF00017_PF00041_5dc4_A1_B1.txt	    25	   150	A   	 169	B   	  73	     12.22	     12.22	      15.6	P00519_150_L	P02751_1517_T	0.11549	618
PF00017_PF00041_4je4_A1_B1.txt	    27	   157	A   	  32	B   	  77	     12.94	      14.2	     19.19	Q06124_32_R	P02751_1524_R	0.10773	711
PF00017_PF00041_5dc4_A1_B1.txt	    27	    80	A   	 171	B   	   6	      21.9	      21.9	      21.9	P00519_152_R	P02751_1450_V	0.10614	724
PF00017_PF00041_5dc4_A1_B1.txt	    27	   146	A   	 171	B   	  69	     26.39	     26.39	     26.39	P00519_152_R	P02751_1513_D	0.10359	769
PF00017_PF00041_5dc0_B1_A1.txt	    24	    97	B   	 168	A   	  21	     17.28	     17.28	     19.89	None	P02751_1467_S	0.10357	770
PF00017_PF00041_5dc9_A1_B1.txt	    74	   110	A   	 218	B   	  33	     23.96	     23.96	     23.96	P00519_199_V	P02751_1477_Y	0.10212	788
PF00017_PF00041_4je4_A1_B1.txt	    27	    97	A   	  32	B   	  21	     28.26	     28.26	     28.26	Q06124_32_R	P02751_1467_S	0.09937	819
PF00017_PF00041_4je4_A1_B1.txt	     5	   139	A   	  10	B   	  61	     22.66	     22.66	     22.66	Q06124_10_N	P02751_1506_S	0.09932	820
PF00017_PF00041_4je4_A1_B1.txt	    27	   140	A   	  32	B   	  62	     24.88	     24.88	     24.88	Q06124_32_R	P02751_1507_G	0.09222	925
PF00017_PF00041_4je4_A1_B1.txt	    11	   105	A   	  16	B   	  26	     27.35	     27.35	     27.35	Q06124_16_A	P02751_1472_A	0.09086	954
```

In some cases, one would like to restrict the number of structures condsidered. This is why the `mind`
argument can also take a path to a text file containing a single column with all desired PDB codes.
As an example, inspect `support_files/template_pdb_list_for_mindist.txt`.
This command will then take into account the first five structures listed above:

`
python3 dca2pdb.py -pf1 PF00017 -pf2 PF00041 -dca ../examples/plmDCA_outputs/PF00017_PF00041_plmdca.txt -mind support_files/template_pdb_list_for_mindist.txt
`

### Testing domain-domain interactions without the structure of the complex

Ultimately, we want to use DCA as a domain-domain contact prediction method, and we are thus interested 
in those cases where the structure of the complex is missing. So far, the only available implementation
is adding the flag "-a" to a simple backmap research:

`
python3 dca2pdb.py -pf PF00940 -pdb 3e2e -a -vmd </path/to/bin/vmd>
`

**Note:** in order to run, this function needs **VMD** Version 1.9.3. In order to link it,
you have to type the path of the VMD executable after the -vmd flag.
The last section of the stdout summary is:

```
SASA:
../examples/results/3e2e_SASA.txt
```

and tells where to find the solvent accessible surface area analysis for the chosen domain.

## Developer's Guide

The DCA-2-PDB suite makes use of the precomputed Pfam interaction benchmark presented in
this repository as well as of several files taken directly from the [Pfam][pfam] database,
which are used in order to faithfully backmap the DCA predictions onto the PDB structures.

The guide to the algorithm is coming soon.

[pfam]: https://pfam.xfam.org/
[vmd]: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD
[plmDCAcode]: https://github.com/pagnani/PlmDCA
[ParalogMatchingcode]: https://github.com/Mirmu/ParalogMatching.jl
