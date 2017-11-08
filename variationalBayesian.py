#! /usr/bin/python
"""

Usage: runVBW.py -k -s layer_lines(simulated)

"""
__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import optparse
import vbwSC

def produce_output(output, file_list):
    #1. Match file list with weight list
    logfile = open(output+".log","w")
    output_lines = open(output).readlines()
    weights = output_lines[-6].split(" ")[:-4]
    output_JSD = output_lines[-5]
    output_chi2 = output_lines[-3]
    output_ModelEvidence = output_lines[-1]

    logfile.write("Log file from variational Bayesian algorithm\n")
    logfile.write(output_chi2)
    logfile.write(output_JSD)
    logfile.write(output_ModelEvidence)

    logfile.write("Models selected:\n")
    structure_file_list = open(file_list).readlines()[0].strip("\n").split(" ")
    if len(weights) != len(structure_file_list):
        raise Exception("Weights and file list have different size!")
    for index, weight in enumerate(weights):
        if float(weight)>0.0:
            logfile.write(structure_file_list[index]+"\n")
    logfile.close()

if __name__=="__main__":
    doc = """
        Python interface to Variational Bayesian algorithm
        Usage: python runVBW.py --help
    """
    print doc
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-r", "--restart", dest="restart",default = 0,
                      type = 'int',
                      help="Restart or not")
    parser.add_option("-n", "--number_of_strcutures", dest="nstruct",default = None,
                      type = 'int',
                      help="Number of strcutural model [OBLIGATORY]")
    parser.add_option("-m", "--measures", dest="measures",default = None,
                      type = 'int',
                      help="Number of measurements on which algorithm is evaluated [OBLIGATORY]")
    parser.add_option("-M", "--cs_measures", dest="cs_measures",default = None,
                      type = 'int',
                      help="Number of chemical shifts measurements  [OBLIGATORY]")
    parser.add_option("-v", "--number_of_curves", dest="ncurves",default = 1,
                      type = 'int',
                      help="Number of scattering curves [OBLIGATORY]")
    parser.add_option("-p", "--priors", dest="priors",
                      help="Prior weights [OBLIGATORY]")
    parser.add_option("-a", "--strcuture_energies", dest="structure_energies", default="None",
                      help="Energies of strcutures used to setup priors")
    parser.add_option("-s", "--simulated", dest="simulated",
                      help="Simulated SAXS curves [OBLIGATORY]")
    parser.add_option("-e", "--experimental", dest="experimental",
                      help="Experimental SAXS curves [OBLIGATORY]")
    parser.add_option("-S", "--cs_simulated", dest="cs_simulated", default="None",
                      help="Simulated CS data [OBLIGATORY]")
    parser.add_option("-E", "--cs_experimental", dest="cs_experimental", default="None",
                      help="Experimental CS data [OBLIGATORY]")
    parser.add_option("-R", "--cs_rms", dest="cs_rms", default="None",
                      help="RMS of simulated CS data [OBLIGATORY]")
    parser.add_option("-o", "--output", dest="output",
                      help="Output file [OBLIGATORY]")
    parser.add_option("-c", "--cores", dest="nprocs",default = 1,
                      type = 'int',
                      help="Number of proccessors")
    parser.add_option("-w", "--weights", dest="weight_cut",default = None,
                      type = 'float',
                      help="Weight cutoff [OBLIGATORY]")
    parser.add_option("-k", "--skip_vbw", dest="skip_vbw",default = 0,
                      type = 'int',
                      help="Skipping VBW step goes to model evidence directly")
    parser.add_option("-f", "--file_list", dest="file_list",
                      help="List of structure files")
    options, args = parser.parse_args()
    vbwSC.run_vbw(options.restart, options.nstruct, options.priors, options.structure_energies,\
                  options.measures, options.cs_measures, options.simulated, options.ncurves, options.experimental,\
                  options.output, options.nprocs, options.weight_cut, options.skip_vbw,\
                  options.cs_simulated, options.cs_rms, options.cs_experimental)

    produce_output(options.output, options.file_list)