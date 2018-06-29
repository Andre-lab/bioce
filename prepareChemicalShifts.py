#! /usr/bin/python
"""
This script reads cs output from SHIFTX2 and combines for BW input

Usage: prepareChemicalShifts.py --help

"""
__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import optparse
import os


class ShiftX2:
	def __init__(self, dir_name, exp_cs_file):
		self.dir_name = dir_name
		self.cs_list = [cs_file for cs_file in os.listdir(dir_name) if cs_file[-3:]==".cs"]
		self.dat_file = "cs.dat"
		self.err_file = open("cs.err","w")
		self.exp_file = open("cs_exp.dat","w")
		self.exp_err = open("cs_exp.err","w")
		self.dat_dict = {}
		self.exp_res_name = {}
		self.exp_shift = {}
		exp_lines = open(exp_cs_file).readlines()
		self.shiftx2_rms = { "N":1.2328, "CA": 0.3836, "CB": 0.5329, "C":0.5096, "H": 0.1711, "HN": 0.2351, "HA": 0.1081,\
		"CD":0.5059, "CD1":1.0113, "CD2":1.812, "CE":0.3690, "CE1":0.9168, "CE2":0.8092, "CG":0.6635, "CG1":0.7835, "CG2":0.9182,\
		"CZ":1.1673, "HA2":0.2732, "HA3":0.3007, "HB":0.1464, "HB2":0.1814, "HB3":0.2109, "HD1":0.1760, "HD2":0.1816, "HD3":0.1796, \
		"HE":0.3342, "HE1":0.1949, "HE2":0.1590, "HE3":0.3219, "HG":0.3826, "HG1":0.1489, "HG12":0.1881, "HG13":0.1269, "HG2":0.1401,\
		"HG3":0.1689, "HZ":0.2610}
		self.exp_errors = {"N":0.3, "H":0.025, "C":0.3}
		for exp_line in exp_lines:
			values = exp_line.split()
			self.exp_shift[(values[0],values[1],values[2])]=values[3]
			if values[0] not in self.exp_res_name:
				self.exp_res_name[values[0]] = [(values[1],values[2])]
			else:
				self.exp_res_name[values[0]].append((values[1],values[2]))


	def read_cs_file(self):
		for count, structure in  enumerate(self.cs_list):
			cs_file_name = structure.strip("\n")
			cs_file = open(os.path.join(self.dir_name,cs_file_name))
			line = 0
			for cs_line in cs_file.readlines():
				if cs_line[:3] == "NUM":
					continue
				values=cs_line.split(",")
				res_number = values[0]
				res_name = values[1]
				atom_name = values[2]
				#check if there is corresponding value in experimental file
				if res_number not in self.exp_res_name.keys():
					continue
				if (res_name,atom_name) not in self.exp_res_name[res_number]:
					continue
				cs =  float(values[3])
				if atom_name not in self.shiftx2_rms.keys():
					#print atom_name
					continue
				if count == 0:
					self.exp_file.write(self.exp_shift[(res_number,res_name,atom_name)]+"\n")
					self.exp_err.write(str(self.exp_errors[atom_name[0]])+"\n")
					err = self.shiftx2_rms[atom_name]
					self.err_file.write(str(err)+"\n")
				if line not in self.dat_dict.keys():
					self.dat_dict[line] = [cs]
				else:
					self.dat_dict[line].append(cs)
				line+=1

	def write_to_files(self):
		data_file = open(self.dat_file, "w")
		for line in self.dat_dict.keys():
			for val in self.dat_dict[line]:
				data_file.write(str(val)+" ")
			data_file.write("\n")

if __name__=="__main__":
	doc = """
	Reads ShiftX2 list of structural files and produces data and error file for BW analysis
	Chemical shifts should be produced with SHIFTX2 before running this script
	Usage: python prepareChemicalShifts.py --help
	"""
	print(doc)
	usage = "usage: %prog [options] args"
	option_parser_class = optparse.OptionParser
	parser = option_parser_class( usage = usage, version='0.1' )

	parser.add_option("-c", "--chemical_shifts_library", dest="dir_name",
                      help="Directory containing pdb files [OBLIGATORY]")
	parser.add_option("-e", "--exp", dest="exp_file",
                      help="Experimental tab file [OBLIGATORY]")

	options, args = parser.parse_args()
	dir_name = options.dir_name
	exp_file = options.exp_file
	shiftx2 = ShiftX2(dir_name, exp_file)
	shiftx2.read_cs_file()
	shiftx2.write_to_files()

