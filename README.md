# PDB_Protein_Analysis

This is a Repo consists of two modules:
1. A PDB (Protein Data Bank) dataset parser that will give cleaned ATOM information.
2. A Ramachandran Analysis tool.

## PDB-Parser
See `PDB_Parser.py`.

It goes to RCSB PDB (Protein Data Bank) and download (stream) the list of non-redundant protein structure files at 30% sequence identity level. 
The resulting text file "clusters-by-entity-30.txt" contains over 300,000 lines, each of which corresponds to a cluster of single-chain sequences and structures (those four alphanumeric characters are PDB IDs and they are followed by "_" and then by a polymer entity identifier, not chain identifier)

Next, loops over the largest 100 clusters (the first 100 lines) in the list, select one random structure for each cluster/line.

Finally, it extracts ATOM information from PDB and FASTA dataset, and returns a cleaned dataframe with: atom_name, residue_name, x, y, z.

## Ramachandran Analysis tool
See `Ramachandran_Analysis.py` and find the experiment result at `Ramachandran_Report.pdf`.

Gives the Ramachandran Plots (scatter plots) for:
 (a) all residues but glycines and prolines
 (b) all glycines
 (c) all prolines

## Experiments
To get the Ramachandran Plots, execute  `$python __main__.py`.

Note that Ramachandran_Analysis must be initialized with a pandas dataframe having the format (columns): atom_name, residue_name, x, y, z.

## Future Work/ Improvements
1. Current version skips the PDBx/mmCIF Format, thus the "first 100" structures is actually giving fewer (94, in the report case).
2. Current version simply accumulate all ATOMs, instead of making them into groups of chains, which may cause mis-calculations (alleviated to 128 occurrences in over 500000 amino acids).
3. Multi-processing/multi-threading methods can be used to improve the speed.

