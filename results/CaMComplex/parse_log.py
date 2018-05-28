import operator
f = open("entire6.txt")
comb_dict = {}
cam_res = 148 
lines = f.readlines()
for line in lines:
    if "Matchmaker" in line:
        comb_name = line
    if "RMSD" in line:
        atom_pairs = float(line.split(" ")[2])
        rmsd = float(line.split(" ")[6])
        score = (cam_res-atom_pairs)*rmsd
        comb_dict[comb_name+line] = score + rmsd
sorted_x = sorted(comb_dict.items(), key=operator.itemgetter(1))
print(sorted_x[:10])
