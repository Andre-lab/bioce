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
import matplotlib
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

def make_intensity_plot(data, simulated_data, log=False):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    matplotlib.rcParams.update({'font.size': 18})

    qvector = data[:,0]
    exp_intensities = data[:,1]
    exp_errors = data[:,2]
    sim_intensities = simulated_data[:,1]

    #plt.plot(qvector, exp_intensities, 'ko', markersize=4, mfc="none")
    plt.errorbar(qvector, exp_intensities, yerr=exp_errors, fmt="ko", markersize=4, alpha=0.4, mfc="none")
    plt.plot(qvector, sim_intensities, 'o', markersize=5)
    plt.yscale('log')

    #plt.legend(handles=[line1])
    #if log:
    #    plt.semilogy(qvector, exp_intensities, 'g-', linewidth=1.5)
    #    plt.semilogy(qvector, sim_intensities, 'g-', linewidth=1.5)
    #else:
    #    plt.plot(qvector, sim_intensities, 'o')
    #    plt.plot(qvector, exp_intensities, 'o')
    #    #plt.plot(qvector, intensities, 'k-', linewidth=1)

    plt.ylabel("$log(Intenisty)$")
    plt.xlabel("$q [\AA^{-1}]$")
    #plt.ylabel("RMSD")
    #plt.xlabel("Noise $(\sigma)$")
    #plt.figure(figsize=(8, 6))
    plt.savefig("SASfit.png", dpi=300, bbox_inches='tight')
    plt.show()


def make_bar_plot():
    """
    Creatiing bar chart
    :param data:
    :param data_name:
    :return:
    """
    matplotlib.rcParams.update({'font.size': 22})
    xvalues = [1,3,5]
    heights = [0.424346, 0.54757, 0.476754]
    widhts = [1.4, 1.4, 1.4]
    fig, ax = plt.subplots(1, 1)
    plt.bar(xvalues, heights, widhts, linewidth=3, edgecolor='orange', color='None')
    ax.set_ylabel("Model Evidence")
    ax.set_xticks(xvalues)
    ax.set_xticklabels(['2models', '3models', '4models'])
    fig.savefig("bar_plot.png", dpi=600, bbox_inches='tight')
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
    exp_errors = data[:,2]
    sim_intensities = combined_data
    line1, = plt.plot(qvector, (exp_intensities - sim_intensities)/exp_errors, 'o')
    #plt.fill_between(qvector, -exp_errors, exp_errors, color='lightyellow')
    #plt.legend(handles=[line1])
    plt.ylabel("$\Delta/\sigma$")
    plt.xlabel("$q[\AA^{-1}]$")
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
    #make_intensity_plot(data, combined_data, data_file)
    make_bar_plot()
