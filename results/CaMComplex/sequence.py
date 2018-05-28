from Bio import SeqIO
handle = open("trans_1CLL_full_noCA_0031_894.pdb", "rU")
for record in SeqIO.parse(handle, "pdb-seqres") :
    print (">" + record.id + "\n" + record.seq)
handle.close()
