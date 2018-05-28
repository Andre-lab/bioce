from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

## The only argument is the PDB file
pdbFile = sys.argv[1]

## First, open and parse the protein file
p = PDBParser(PERMISSIVE=1)
structure = p.get_structure(pdbFile, pdbFile)

## Now go through the hierarchy of the PDB file
##
## 1- Structure
##      2- Model
##          3- Chains
##              4- Residues
##

for model in structure:
    for chain in model:
        seq = list()
        chainID = chain.get_id()

        for residue in chain:
            ## The test below checks if the amino acid
            ## is one of the 20 standard amino acids
            ## Some proteins have "UNK" or "XXX", or other symbols
            ## for missing or unknown residues
            if is_aa(residue.get_resname(), standard=True):
                seq.append(three_to_one(residue.get_resname()))
            else:
                seq.append("X")

        ## This line is used to display the sequence from each chain

        print ">Chain_" + chainID + "\n" + str("".join(seq))

        ## The next two lines create a sequence object/record
        ## It might be useful if you want to do something with the sequence later

        myProt = Seq(str(''.join(seq)), IUPAC.protein)
        seqObj = SeqRecord(myProt, id=chainID, name="", description="")

## The end