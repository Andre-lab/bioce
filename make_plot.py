"""
Full bayesian inference using stan
"""
from __future__ import print_function

__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import os
import optparse

import numpy as np
import matplotlib.pyplot as plt

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

def make_plot(data, data_name, log=False):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    qvector = data[:,0]
    #intensities = data[:,1:]

    intensities = data[:,1]
    line1, = plt.plot(qvector, intensities, '-o')
    plt.legend(handles=[line1])
    #if log:
    #    plt.semilogy(qvector, intensities, 'g-', linewidth=1.5)
    #else:
    #    plt.plot(qvector, intensities, '-o')
    #    #plt.plot(qvector, intensities, 'k-', linewidth=1)

    #plt.ylabel("$log(Intenisty)$")
    #plt.xlabel("$q [\AA]$")
    plt.ylabel("RMSD")
    plt.xlabel("Noise $(\sigma)$")
    plt.savefig(data_name+".png", dpi=600)
    plt.show()

def make_bar_plot(data, data_name):
    """
    Creatiing bar chart
    :param data:
    :param data_name:
    :return:
    """
    xvalues = [1,3,5]
    heights = [0.424346, 0.54757, 0.476754]
    widhts = [1.8, 1.8, 1.8]
    fig, ax = plt.subplots(1, 1)
    plt.bar(xvalues, heights, widhts)
    ax.set_ylabel("Model Evidence")
    ax.set_xticks(xvalues)
    ax.set_xticklabels(['2models', '3models', '4models'])
    fig.savefig(data_name+".png", dpi=600)
    plt.show()
def make_residue_plot(data, combined_data, data_name, log=False):
    """
    Making residual plot. Difference between experimental and simulated curves
    Scaling factor needs to be applied first
    Can be compared with teh experimental error
    :param data:
    :param log:
    :return:
    """
    qvector = data[:,0]
    exp_intensities = data[:,1]
    sim_intensities = combined_data
    line1, = plt.plot(qvector, exp_intensities - sim_intensities, 'o')
    #plt.legend(handles=[line1])
    plt.ylabel("Residuals")
    plt.xlabel("$q[\AA]$")
    plt.savefig(data_name+".png", dpi=600)
    plt.show()
if __name__=="__main__":
    doc = """
        Python interface to produce sacttering plots.
        The script reads X vector and Y matrix containing one or more intensity set
        The available options should include scattering plots, connected or
        disconnected, on the log or normal scale
    """
    print(doc)
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-d", "--data", dest="data_file", default=None,
                      help="Data file containing qvec, intensities [OBLIGATORY]")
    parser.add_option("-c", "--cmb_data", dest="combined_data_file", default=None,
                      help="Data file containing qvec, intensities [OBLIGATORY]")
    options, args = parser.parse_args()

    data_file = options.data_file
    combined_data_file = options.combined_data_file
    data = read_file_safe(data_file)
    combined_data = read_file_safe(combined_data_file)
    make_residue_plot(data, combined_data, data_file)

