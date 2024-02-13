# __main__.py

from imports import *

# Import custom classes
from classes.class1 import PDB_Parser
from classes.class2 import Ramachandran_Analysis

def main():
    # new a PDB_Parser and get the atom data from PDB/FASTA DB
    p = PDB_Parser()
    atom_data = p.main()
    
    # new a Ramachandran_Analysis and plot the Ramachandran Plot for three subsets
    r = Ramachandran_Analysis(atom_data)
    r.main()

if __name__ == "__main__":
    main()