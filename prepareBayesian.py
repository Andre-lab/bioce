#! /usr/bin/python
"""

Usage: runVBW.py -k -s layer_lines(simulated)

"""
__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import os
import shlex
import os.path
import string
from subprocess import Popen, PIPE
import optparse

def generate_scattering_profiles(dir_name, exp_file):
	file_list = os.listdir(dir_name)
	out_file = open("SimulatedIntensities.txt","w")
	weights = open("weights.txt","w")
	structure_list = open("strcutures.txt","w")
	int_dict = {}
	name_dict = {}

	for count, pdb in enumerate(file_list):
		pdb = pdb.strip("\n")
		if pdb[-4:] != ".pdb":
			continue
		cmd_line = "foxs "+os.path.join(dir_name,pdb)+" "+exp_file
		proc = Popen(shlex.split(cmd_line), stdout=PIPE, shell=False)
		os.waitpid(proc.pid,0)
		out, err = proc.communicate()
		foxs_file = dir_name+"/"+pdb[:-4]+"_"+exp_file
		int_dict[count] = []
		name_dict[count] = pdb[:-4]
		foxs_lines = open(foxs_file).readlines()[1:]
		norm_c = 1E8
		#TODO: Check this scaling parameter
		scaling_c = norm_c*float(foxs_lines[0].split(",")[1][12:])
		#print("C scaling", scaling_c)
		for line in foxs_lines[2:]:
			words=line.split()
			intensity_foxs = float(words[2])/scaling_c
			int_dict[count].append(str(intensity_foxs))

	number_of_structures = len(int_dict.keys())

	for count in int_dict.keys():
		weights.write(str(1.0/number_of_structures)+" ")
		structure_list.write(name_dict[count]+" ")
	structure_list.write("\n")
	weights.write("\n")

	for line in range(len(foxs_lines[2:])):
		for count in int_dict.keys():
			out_file.write(int_dict[count][line]+" ")
		out_file.write("\n")

	structure_list.close()
	weights.close()
	out_file.close()

if __name__=="__main__":
	doc = """
        Script to generate scattering curves from pdb files
        "FOXS" program must be installed in PATH
        Usage: python prepareBayesian.py --help
    """
	print(doc)
	usage = "usage: %prog [options] args"

	if not any([os.path.exists(os.path.join(p, "foxs"))
							   for p in os.environ["PATH"].split(os.pathsep)]):
		raise RuntimeError("FoXS not avaialble in the PATH."
						   "Please add FoXS locationn to PATH and run script again")
	option_parser_class = optparse.OptionParser
	parser = option_parser_class( usage = usage, version='0.1' )
	parser.add_option("-s", "--structure_library", dest="dir_name",
                      help="Directory containing pdb files")
	parser.add_option("-e", "--experimental", dest="exp_file",
                      help="Experimental SAXS curves [OBLIGATORY]")
	options, args = parser.parse_args()
	generate_scattering_profiles(options.dir_name,options.exp_file)