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
import pystan
import psisloo
import matplotlib.pyplot as plt

from statistics import calculateChiCrysol, calculateChemShiftsChi, JensenShannonDiv
from stan_models import stan_code, stan_code_CS, stan_code_EP, stan_code_EP_CS, \
    psisloo_quanities

def execute_stan(experimental, simulated, priors, iterations, chains, njobs):
    """

    :param experimental:
    :param simulated:
    :param priors:
    :param iterations:
    :param chains:
    :param njobs:
    :return:
    """
    print("Starting SAXS inference")
    stan_dat = {"sim_curves": simulated,
            "target_curve": experimental[:,1],
            "target_errors": experimental[:,2],
            "n_measures" : np.shape(experimental)[0],
            "n_structures" : np.shape(simulated)[1],
            "priors":priors}
    sm = pystan.StanModel(model_code=stan_code)
    fit = sm.sampling(data=stan_dat, iter=iterations, chains=chains, n_jobs=njobs, sample_file="saved_samples.txt")

    fig = fit.plot(pars="weights")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_weights.png", dpi=300)

    fig = fit.plot(pars="scale")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_scale.png", dpi=300)

    return fit


def execute_stan_EP(experimental, simulated, priors, iterations, chains, njobs):
    """

    :param experimental:
    :param simulated:
    :param priors:
    :param iterations:
    :param chains:
    :param njobs:
    :return:
    """
    print("Starting SAXS+EP inference")
    stan_dat = {"sim_curves": simulated,
            "target_curve": experimental[:,1],
            "target_errors": experimental[:,2],
            "n_measures" : np.shape(experimental)[0],
            "n_structures" : np.shape(simulated)[1],
            "energy_priors":priors}
    sm = pystan.StanModel(model_code=stan_code_EP)
    fit = sm.sampling(data=stan_dat, iter=iterations, chains=chains, n_jobs=njobs, sample_file="saved_samples.txt")

    fig = fit.plot(pars="weights")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_weights_EP.png", dpi=300)

    fig = fit.plot(pars="boltzman_shift")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_boltzmann_EP.png", dpi=300)

    fig = fit.plot(pars="scale")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_scale_EP.png", dpi=300)


    return fit

def execute_stan_EP_CS(experimental, simulated, priors,
                       simulated_cs, simulated_cserr, experimental_cs,
                       iterations, chains, njobs):
    """

    :param experimental:
    :param simulated:
    :param priors:
    :param simulated_cs:
    :param simulated_cserr:
    :param experimental_cs:
    :param iterations:
    :param chains:
    :param njobs:
    :return:
    """
    print("Starting SAXS+CS+EP inference")
    stan_dat = {"sim_saxs": simulated,
            "target_saxs": experimental[:,1],
            "target_saxserr": experimental[:,2],
            "sim_css": simulated_cs,
            "sim_cserr": simulated_cserr,
            "target_cs": experimental_cs[:,0],
            "target_cserr": experimental_cs[:,1],
            "n_measures" : np.shape(experimental)[0],
            "m_measures" : np.shape(experimental_cs)[0],
            "n_structures" : np.shape(simulated)[1],
            "energy_priors":priors}

    sm = pystan.StanModel(model_code=stan_code_EP_CS)
    fit = sm.sampling(data=stan_dat, iter=iterations, chains=chains, n_jobs=njobs, sample_file="saved_samples.txt")

    fig = fit.plot(pars="weights")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_weights_EP_CS.png", dpi=300)

    fig = fit.plot(pars="boltzman_shift")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_boltzmann_EP_CS.png", dpi=300)

    fig = fit.plot(pars="scale")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_scale_EP_CS.png", dpi=300)

    return fit

def execute_stan_CS(experimental, simulated, priors,
                       simulated_cs, simulated_cserr, experimental_cs,
                       iterations, chains, njobs):
    """

    :param experimental:
    :param simulated:
    :param priors:
    :param simulated_cs:
    :param simulated_cserr:
    :param experimental_cs:
    :param iterations:
    :param chains:
    :param njobs:
    :return:
    """
    print("Starting SAXS+CS inference")
    stan_dat = {"sim_saxs": simulated,
            "target_saxs": experimental[:,1],
            "target_saxserr": experimental[:,2],
            "sim_css": simulated_cs,
            "sim_cserr": simulated_cserr,
            "target_cs": experimental_cs[:,0],
            "target_cserr": experimental_cs[:,1],
            "n_measures" : np.shape(experimental)[0],
            "m_measures" : np.shape(experimental_cs)[0],
            "n_structures" : np.shape(simulated)[1],
            "priors":priors}

    sm = pystan.StanModel(model_code=stan_code_CS)
    fit = sm.sampling(data=stan_dat, iter=iterations, chains=chains, n_jobs=njobs, sample_file="saved_samples.txt")

    fig = fit.plot(pars="weights")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_weights_CS.png", dpi=300)

    fig = fit.plot(pars="scale")
    fig.subplots_adjust(wspace=0.8)
    fig.savefig("stan_scale.png", dpi=300)
    return fit

def execute_bws(experimental, simulated, priors, file_names, threshold,
                iterations, chains, njobs):
    """
    Bayesian weighting with selection algorithm
    Weight threshold is taken as an input. Several simulations are
    run during each iteration as model comparison doesn't always decrease

    :param experimental:
    :param simulated:
    :param priors:
    :param names_file:
    :param iterations:
    :param chains:
    :param njobs:
    :return:
    """

    sim_curves = simulated
    target_curve = experimental[:,1]
    target_errors = experimental[:,2]
    n_measures = np.shape(experimental)[0]
    n_structures = np.shape(simulated)[1]
    alphas = priors
    log_file = open("full_bayesian.log","w")

    post_loo = None
    last_loo = None
    model_comp_diff = 1
    sm = pystan.StanModel(model_code=stan_code+psisloo_quanities, iter=iterations, chains=chains,
                          n_jobs=njobs)

    iteration = 0
    repeat_iteration = 0
    while (model_comp_diff > 0 or repeat_iteration < 5):
        post_loo = last_loo
        log_file.write("Starting iteration "+str(iteration)+" with "
                   +str(n_structures)+" models\n")
        stan_dat = {"sim_curves": sim_curves,
            "target_curve": target_curve,
            "target_errors": target_errors,
            "n_measures" : n_measures,
            "n_structures" : n_structures,
            "priors" : alphas}

        fit = sm.sampling(data=stan_dat, iter=iterations, chains=chains,
                          n_jobs=njobs)

        #Calculating psis loo
        stan_chain=fit.extract()
        current_loo = psisloo.psisloo(stan_chain['loglikes'])
        log_file.write("greater than 0.5 ")
        log_file.write(str(current_loo.print_summary()[0])+"\n")
        log_file.write("greater than 1 ")
        log_file.write(str(current_loo.print_summary()[1])+"\n")
        if last_loo:
            log_file.write("Model comparison: \n")
            model_comp = psisloo.loo_compare(last_loo,current_loo)
            model_comp_diff = model_comp['diff']
            log_file.write(str(model_comp['diff'])+"\n")
            log_file.write(str(model_comp['se_diff'])+"\n")

        #TODO check how many time condition has been met and exit if it doesn't improve over 5
        if last_loo and model_comp_diff < 0:
            repeat_iteration += 1
            log_file.write("\nModel comp < 0 repeating iteration\n")
            log_file.write("Repeating "+str(repeat_iteration)+" time\n")
            continue

        #Storing data for next simulation
        last_loo = current_loo
        repeat_iteration = 0
        current_weights = fit.summary()['summary'][:,0][:n_structures]
        sim_curves = sim_curves[:,current_weights>threshold]
        alphas = alphas[current_weights>threshold]
        n_structures = np.shape(sim_curves)[1]
        file_names = file_names[current_weights>threshold]
        #np.savetxt("fit.txt",fit.summary()['summary'][:,0],delimiter=" ")
        fig = fit.plot()
        fig.subplots_adjust(wspace=0.8)
        fig.savefig("stan_fit_"+str(iteration)+".png",  dpi=300)
        bayesian_weights, jsd, crysol_chi2 = calculate_stats(fit,
                                        experimental, simulated)
        log_file.write("JSD: "+str(jsd))
        log_file.write("\nchi2: "+str(crysol_chi2))
        iteration += 1

    log_file.write("Final Model comparison: \n")
    model_comp = psisloo.loo_compare(post_loo,last_loo)
    log_file.write(str(model_comp['diff'])+"\n")
    log_file.write(str(model_comp['se_diff'])+"\n")

    bayesian_weights, jsd, crysol_chi2 = calculate_stats(fit,
                                        experimental, sim_curves)
    log_file.write("Final weights:\n")
    for index, weight in enumerate(bayesian_weights):
        log_file.write(str(file_names[index])+" : "+str(weight)+"\n")

    log_file.write("\nJSD: "+str(jsd)+"\n")
    log_file.write("\nchi2: "+str(crysol_chi2)+"\n")
    log_file.close()

# def execute_bws_EP_CS(experimental, simulated, priors,
#                        simulated_cs, simulated_cserr, experimental_cs,
#                        iterations, chains, njobs):
#     """
#
#     :param experimental:
#     :param simulated:
#     :param priors:
#     :param simulated_cs:
#     :param simulated_cserr:
#     :param experimental_cs:
#     :param iterations:
#     :param chains:
#     :param njobs:
#     :return:
#     """

def calculate_stats(fit, experimental, simulated, cs_simulated=None,
                    cs_rms=None, cs_experimental=None):
    """
    Calculates statistics based on stan model
    :param fit:
    :return:
    """
    #la = fit.extract(permuted=True)  # return a dictionary of arrays
    #mu = la['weights']

    ## return an array of three dimensions: iterations, chains, parameters
    results_array = fit.extract(permuted=False, inc_warmup=False)

    nsamples = 0
    jsd_sum = 0.0
    bayesian_weights = np.zeros(np.shape(simulated)[1])
    for iteration in results_array:
        for parameters in iteration:
            current_weights = parameters[:np.shape(simulated)[1]]
            bayesian_weights+=current_weights
            nsamples+=1
    bayesian_weights=bayesian_weights/nsamples

    for iteration in results_array:
        for parameters in iteration:
            current_weights = parameters[:np.shape(simulated)[1]]
            jsd_sum+=JensenShannonDiv(current_weights, bayesian_weights)
    jsd = (np.sqrt(jsd_sum/nsamples))


    crysol_chi2 = calculateChiCrysol(np.dot(bayesian_weights,
                            np.transpose(simulated)), experimental[:,1],
                            experimental[:,2])
    try:
        if cs_experimental.any() != None:
            chemical_shifts_on = True
    except:
        chemical_shifts_on = False

    if chemical_shifts_on:
        chemshift_chi2 = calculateChemShiftsChi(np.dot(bayesian_weights,
                            np.transpose(cs_simulated)), cs_experimental[:,0],
                            cs_experimental[:,1], cs_rms)
        return bayesian_weights, jsd, crysol_chi2, chemshift_chi2
    else:
        return bayesian_weights, jsd, crysol_chi2

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

if __name__=="__main__":
    doc = """
        Python interface to Variational Bayesian algorithm
        Usage: python runVBW.py --help
    """
    print(doc)
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-p", "--priors", dest="priors_file", default=None,
                      help="Prior weights [OBLIGATORY]")
    parser.add_option("-P", "--strcuture_energies", dest="structure_energies", default=None,
                      help="Energies of strcutures used to setup priors")
    parser.add_option("-s", "--simulated", dest="simulated_file",
                      help="Simulated SAXS curves [OBLIGATORY]")
    parser.add_option("-e", "--experimental", dest="experimental_file",
                      help="Experimental SAXS curves [OBLIGATORY]")
    parser.add_option("-S", "--cs_simulated", dest="cs_simulated_file", default=None,
                      help="Simulated CS data [OBLIGATORY]")
    parser.add_option("-E", "--cs_experimental", dest="cs_experimental_file", default=None,
                      help="Experimental CS data [OBLIGATORY]")
    parser.add_option("-R", "--cs_rms", dest="cs_rms_file", default=None,
                      help="RMS of simulated CS data [OBLIGATORY]")
    parser.add_option("-o", "--output", dest="output",
                      help="Output file [OBLIGATORY]")
    parser.add_option("-w", "--weights", dest="weight_cut",default = None,
                      type = 'float',
                      help="Weight cutoff [OBLIGATORY]")
    parser.add_option("-f", "--file_names", dest="names_file",default = None,
                      help="File containing names of structures to process")
    parser.add_option("-i", "--iterations", dest="iterations",default = 2000,
                      type = 'int',
                      help="Number of iterations")
    parser.add_option("-j", "--jobs", dest="njobs",default = 8,
                      type = 'int',
                      help="Number of procceses to run")
    parser.add_option("-c", "--chains", dest="chains",default = 4,
                      type = 'int',
                      help="Number of chains")
    options, args = parser.parse_args()

    njobs = options.njobs
    iterations = options.iterations
    chains = options.chains

    #Files loading
    experimental = read_file_safe(options.experimental_file)
    simulated =read_file_safe(options.simulated_file)

    #Loading file names
    if options.names_file:
        file_names = read_file_safe(options.names_file,'unicode')

    if options.priors_file:
        priors = read_file_safe(options.priors_file)
    else:
        priors = read_file_safe(options.structure_energies)

    if options.cs_simulated_file:
        cs_simulated = read_file_safe(options.cs_simulated_file)
        cs_experimental = read_file_safe(options.cs_experimental_file)
        cs_rms = read_file_safe(options.cs_rms_file)

    #If weight threshold set running Bayesian with selection
    if options.weight_cut:
        execute_bws(experimental, simulated, priors, file_names, options.weight_cut,
                iterations, chains, njobs)
    else:
        #Running siumlations depending on options selected
        if options.structure_energies:
            if options.cs_simulated_file:
                fit = execute_stan_EP_CS(experimental, simulated, priors,
                                         cs_simulated, cs_rms, cs_experimental,
                                         iterations, chains, njobs)
            else:
                fit = execute_stan_EP(experimental, simulated, priors,
                                      iterations, chains, njobs)
        else:
            if options.cs_simulated_file:
                fit = execute_stan_CS(experimental, simulated, priors,
                                      cs_simulated, cs_rms, cs_experimental,
                                      iterations, chains, njobs)
            else:
                fit = execute_stan(experimental, simulated, priors,
                                   iterations, chains, njobs)

        if options.cs_simulated_file:
            bayesian_weights, jsd, crysol_chi2, chemshift_chi2 = \
                calculate_stats(fit, experimental, simulated,
                                cs_simulated, cs_rms, cs_experimental)
            for index,fname in enumerate(file_names):
                print(fname,bayesian_weights[index])
            print("JSD: "+str(jsd))
            print("Chi2 SAXS:"+str(crysol_chi2))
            print("Chi2 CS:"+str(chemshift_chi2))
        else:
            bayesian_weights, jsd, crysol_chi2 = \
                calculate_stats(fit, experimental, simulated)
            for index,fname in enumerate(file_names):
                print(fname,bayesian_weights[index])
            print("JSD: "+str(jsd))
            print("Chi2 SAXS:"+str(crysol_chi2))
        print(fit)
