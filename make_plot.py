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
        results = np.genfromtxt(filename, dtype=dtype, delimiter=";")
    except IOError as err:
        print(os.strerror(err.errno))
    return results

def make_fig1c(data, data_name, log=False):
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
    #plt.legend(handles=[line1])
    #if log:
    #    plt.semilogy(qvector, intensities, 'g-', linewidth=1.5)
    #else:
    #    plt.plot(qvector, intensities, '-o')
    #    #plt.plot(qvector, intensities, 'k-', linewidth=1)

    #plt.ylabel("$log(Intenisty)$")
    #plt.xlabel("$q [\AA]$")
    plt.ylabel("RMSD")
    plt.xlabel("Noise $(\sigma)$")
    plt.xticks([1,2,3,4,5])
    plt.savefig(data_name+".png", dpi=600, bbox_inches='tight')
    plt.show()


def make_fig1b(data, data_name, log=False):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    matplotlib.rcParams.update({'font.size': 14})
    qvector = data[:,0]
    #intensities = data[:,1:]
    fig = plt.figure()
    ax = plt.subplot(111)

    intensities = data[:,1]
    line1, = ax.plot(qvector, intensities, '-o', label='Ng')
    intensities = data[:,2]
    line2, = ax.plot(qvector, intensities, '-o', label='models recovered')
    intensities = data[:,3]
    line3, = ax.plot(qvector, intensities, '-o', label='models recovered from preset')
    #plt.legend(handles=[line1, line2, line3])

    ax.legend(loc='upper center', bbox_to_anchor=(0.7, 1.17),
              ncol=1, fancybox=True, shadow=True)
    #if log:
    #    plt.semilogy(qvector, intensities, 'g-', linewidth=1.5)
    #else:
    #    plt.plot(qvector, intensities, '-o')
    #    #plt.plot(qvector, intensities, 'k-', linewidth=1)

    #plt.ylabel("$log(Intenisty)$")
    #plt.xlabel("$q [\AA]$")
    plt.ylabel("Ng/Number of models")
    plt.xlabel("Noise $(\sigma)$")
    plt.savefig(data_name+".png", dpi=600)
    plt.show()

def make_s1fig(data, data_name):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    xvalues = data[:,0]

    #SuppFig 1
    for i in range(5):
        yvalues = data[:,i+1]
        line1, = plt.plot(xvalues, yvalues, '-o',  markersize=2)

    plt.ylabel("weights")
    plt.xlabel("Iteration")
    plt.savefig(data_name+".png", dpi=600)
    plt.show()

def make_plot(data, data_name):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    xvalues = data[:,0]

    #SuppFig 1
    #for i in range(5):
    #    yvalues = data[:,i+1]
    #    line1, = plt.plot(xvalues, yvalues, '-o',  markersize=2)

    #SuppFig2
    yvalues = data[:,1]
    line1, = plt.plot(xvalues, yvalues, '-o',  markersize=4, label="energy (on)")
    yvalues = data[:,2]
    line2, = plt.plot(xvalues, yvalues, '-o',  markersize=4, label="energy (off)")

    first_legend = plt.legend(handles=[line1], loc=1)
    # Add the legend manually to the current Axes.
    ax = plt.gca().add_artist(first_legend)

    # Create another legend for the second line.
    plt.legend(handles=[line2], loc=4)
    yint=[2,3,4,5,6]
    plt.yticks(yint)
    plt.ylabel("Models recovered from preset")
    plt.xlabel("Noise $(\sigma$)")
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

    plt.plot(qvector, exp_intensities, 'ko', markersize=4, mfc="none")
    plt.plot(qvector, sim_intensities, 'o', markersize=4, zorder=5)
    plt.errorbar(qvector, exp_intensities, yerr=exp_errors,
                 fmt="ko", markersize=6, mfc='none', alpha=0.6, zorder=0)

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


def make_hbv_plot(data_sets, log=False):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    matplotlib.rcParams.update({'font.size': 18})

    legend_lines = []
    data_labels = ['2mer','4mer','6mer','10mer','full']
    for i, data_file in enumerate(data_sets):
        data = np.genfromtxt(data_file)
        print(np.shape(data))
        qvector = data[:,0]
        exp_intensities = data[:,1]
        line, = plt.plot(qvector, exp_intensities, label=data_labels[i])
        legend_lines.append(line)
    plt.legend(handles=legend_lines)
    plt.yscale('log')
    plt.ylabel("$log(Intenisty)$")
    plt.xlabel("$q [\AA^{-1}]$")
    #plt.ylabel("RMSD")
    #plt.xlabel("Noise $(\sigma)$")
    #plt.figure(figsize=(8, 6))
    plt.savefig("HBV_profiles.png", dpi=300, bbox_inches='tight')
    plt.show()


def make_ppc_plot(data, log=False):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    matplotlib.rcParams.update({'font.size': 18})

    qvector = data[:,0]
    exp_intensities = data[:,1]
    #exp_errors = data[:,3]
    plt.plot(qvector, exp_intensities, 'o', markersize=2)
    sim_intensities = data[:,5]
    plt.plot(qvector, sim_intensities, 'o', markersize=2)
    sim_intensities = data[:,7]
    plt.plot(qvector, sim_intensities, 'o', markersize=2)
    sim_intensities = data[:,9]
    plt.plot(qvector, sim_intensities, 'o', markersize=2)
    sim_intensities = data[:,11]
    plt.plot(qvector, sim_intensities, 'o', markersize=2)
    #plt.ylim(2e-4,4e-4)
    #plt.errorbar(qvector, exp_intensities, yerr=exp_errors,
    #             fmt="ko", markersize=2, mfc='none', alpha=0.6, zorder=0)

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
    widhts = [1.2, 1.2, 1.2]
    fig, ax = plt.subplots(1, 1)
    barlist = plt.bar(xvalues, heights, widhts, linewidth=3, edgecolor='black')
    barlist[0].set_color('blue')
    barlist[0].set_edgecolor('black')
    #barlist[0].set_hatch('/')
    barlist[1].set_color('orange')
    barlist[1].set_edgecolor('black')
    #barlist[1].set_hatch('/')
    barlist[2].set_color('green')
    barlist[2].set_edgecolor('black')
    #barlist[2].set_hatch('/')
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
    data = read_file_safe(data_file)
    data_name = 'fig1b'
    #combined_data_file = options.combined_data_file
    #combined_data = read_file_safe(combined_data_file)
    #make_ppc_plot(data,data_file)
    make_fig1b(data, data_name)
    #data_sets = ['1qgt_2mer00.int','1qgt_4mer00.int','1qgt_6mer00.int','1qgt_10mer00.int','1qgt_full00.int']
    #data_sets = ['1qgt_2mer.pdb.dat','1qgt_4mer.pdb.dat','1qgt_6mer.pdb.dat','1qgt_10mer.pdb.dat','1qgt_full.pdb.dat']
    #make_hbv_plot(data_sets, True)