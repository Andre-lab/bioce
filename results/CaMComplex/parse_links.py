import os
f = open("links.txt")
pdb_ids = []
for line in f.readlines():
    if "http://www.rcsb.org/pdb/search/structidSearch.do?structureId=" in line:
        index = line.find("http://www.rcsb.org/pdb/search/structidSearch.do?structureId=")
        pdb_id = line[index+61:index+65]
        command = "wget http://www.rcsb.org/pdb/files/"+pdb_id+".pdb.gz"
        os.system(command)
        print(pdb_id)
