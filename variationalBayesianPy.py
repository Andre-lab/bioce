#! /usr/bin/python
"""

Usage: runVBW.py -k -s layer_lines(simulated)

"""
__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import numpy as np
from scipy.special import gamma, digamma
from scipy import optimize
import optparse
#import vbwSC
import os


def read_file_safe(filename, dtype="float64"):
    """
    Simple check if file exists
    :param filename:
    :return:
    """
    try:
        results = np.genfromtxt(filename, dtype=dtype)
    except IOError as err:
        print(os.strerror(err.errno))
    return results

def produce_final_output(output, file_list):
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

def read_file_safe(filename, skip_lines, dtype="float64"):
    """
    Simple check if file exists
    :param filename:
    :return:
    """
    try:
        results = np.genfromtxt(filename, dtype=dtype, skip_header=skip_lines)
    except IOError as err:
        print(os.strerror(err.errno))
    return results

def extract_parameters(filename):
    data = read_file_safe(filename,0)
    return np.shape(data)

def setup_constants(simulated, errors):
    """

    :return:
    """
    (number_of_measurements, number_of_structures) = extract_parameters(simulated)

    gamma_over_lhalf = 1.0/gamma(0.5*number_of_measurements)
    gamma_half = gamma(0.5)

    #Matrix with i,j,k indexes
    mixed_terms = np.zeros(number_of_structures,
                           number_of_structures,
                           number_of_measurements)
    #TODO: Do smarter summation here to avoid iterations over the loop
    #TODO: Errors need to have proper shape (reshape may be used).
    for i in range(number_of_structures):
        for j in range(number_of_structures):
            for k in range(number_of_measurements):
                mixed_terms[i,j,k] = simulated[k,i]*simulated[k,j]\
                                     /(errors[k]**errors[k])
    return gamma_over_lhalf, gamma_half, mixed_terms

def Lfunction(alphas, *params):
    """

    :return:
    """
    experimental, errors, intensities, gamma_over_lhalf, gamma_half, mixed_terms = params
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0
    alpha_zero = np.sum(alphas)
    for alpha in alphas:
        term1 += np.log(gamma_half/gamma(alpha))
        #TODO: This lime can be optimizied
        term1 += (alpha - 0.5)*(digamma(alpha)-digamma(alpha_zero))
    for i in range(intensities):
        for j, alpha in enumerate(alphas):
            term2 += np.pow(experimental - intensities[i,j]*alpha/alpha_zero)
    term2 = 0.5*term2/errors

    for i in range(alphas):
        for j in range(alphas):
            alpha_i = alphas[i]
            alpha_j = alphas[j]
            for k in range(intensities):
                if i == j:
                    term3 += mixed_terms[i,j,k] * (alpha_i * (alpha_zero - alpha_i) - alpha_i * alpha_j)
                else:
                    term3 += mixed_terms[i, j, k] * (alpha_zero**2*(alpha_zero + 1))
    term3 = 0.5 * term3
    Lfunc = np.log(gamma_over_lhalf*gamma(alpha_zero)) + term1 + term2 + term3
    return Lfunc

def simulated_annealing(simulated, priors, intensities, errors, gamma_over_lhalf, gamma_half, mixed_terms):
    """
    Performs simulated anealing and returns weights. It is paramterized with alphas
    :return:
    """

    params = intensities, errors, simulated, gamma_over_lhalf, gamma_half, mixed_terms
    alphas0 = priors
    res = optimize.anneal(Lfunction, alphas0, args=params, schedule='boltzmann',
                          full_output=True, maxiter=500, lower=-10,
                          upper=10, dwell=250, disp=True)

def run_variational(simulated_file, priors_file, experimental_file, output_file, file_list, weight_cut):
    """

    :param simulated_file:
    :param priors_file:
    :param experimental_file:
    :param output_file:
    :param file_list:
    :param weight_cut:
    :return:
    """


    experimental = read_file_safe(experimental_file)
    intensities = experimental[:, 1]
    errors = experimental[:,2]
    simulated = read_file_safe(simulated_file)
    priors = read_file_safe(priors_file)

    gamma_over_lhalf, gamma_half, mixed_terms = setup_constants(simulated, errors)

    #TODO: Write to the same format or may them talk to each other
    #TODO: Here we need a loop for model selection
    simulated_annealing(simulated, priors, intensities, errors,
                        gamma_over_lhalf, gamma_half, mixed_terms)

    produce_final_output(output_file, file_list)


if __name__=="__main__":
    doc = """
        Python interface to Variational Bayesian algorithm
        Usage: python variationalBayesian --help
    """
    print(doc)
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-s", "--simulated", dest="simulated",
                      help="Simulated SAXS curves [OBLIGATORY]")
    parser.add_option("-e", "--experimental", dest="experimental",
                      help="Experimental SAXS curves [OBLIGATORY]")
    parser.add_option("-p", "--priors", dest="priors",
                      help="Prior weights [OBLIGATORY]")
    parser.add_option("-o", "--output", dest="output",
                      help="Output file [OBLIGATORY]")
    parser.add_option("-w", "--weights", dest="weight_cut",default = None,
                      type = 'float',
                      help="Weight cutoff [OBLIGATORY]")
    parser.add_option("-f", "--file_list", dest="file_list",
                      help="List of structure files")
    options, args = parser.parse_args()

    run_variational(options.simulated, options.priors, options.experimental,
                    options.output, options.file_list, options.weight_cut)