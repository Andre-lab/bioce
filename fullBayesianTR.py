"""
Full bayesian inference using stan
"""
from __future__ import print_function

__author__ = "Wojtek Potrzebowski"
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import os
import optparse
import pickle

import numpy as np
import pystan
import psisloo
import matplotlib.pyplot as plt
import stan_utility

from statistics import (
    calculateChiCrysol,
    calculateChemShiftsChi,
    JensenShannonDiv,
    calculateScale,
    core_shell_model,
)
from stan_models import stan_code


def combine_curve(experimental, simulated, weights, iteration, fit_name):
    """

    :param simulated:
    :param weights:
    :param scale:
    :return:
    """
    weights = np.array(weights)
    q_column = experimental[:, 0]
    exp_intensities = experimental[:, 1]
    exp_errors = experimental[:, 2]
    weighted_intensities = np.dot(weights, np.transpose(simulated))
    scale = calculateScale(weighted_intensities, exp_intensities, exp_errors)
    combined = scale * weighted_intensities
    np.savetxt(
        os.path.join("iteration_" + str(iteration), fit_name),
        np.transpose((q_column, exp_intensities, combined, exp_errors)),
    )
    # Differential curve
    np.savetxt(
        os.path.join("iteration_" + str(iteration), "diff_" + fit_name),
        np.transpose((q_column, exp_intensities - combined, exp_errors)),
    )


def average_basis_spectra(
    exp_file_list, simulated_intensities, number_of_states, iter, magdiv
):
    """
    Taking all difference weights and corresponding weights or magnitude to obtain basis spectra
    :return:
    """
    last_weight_index = number_of_states + 1
    # Read files from previous iteration
    out_file = read_file_safe(os.path.join("iteration_" + str(iter - 1), "batch_test"))
    combined_intensity = []
    combined_errors = []
    for index, exp_file in enumerate(exp_file_list):
        diff_fit_file = "diff_" + exp_file[:-4] + "fit"
        diff_file = read_file_safe(
            os.path.join("iteration_" + str(iter - 1), diff_fit_file)
        )
        qvector = diff_file[:, 0]
        diff_simulated = diff_file[:, 1]
        diff_errors = diff_file[:, 2]

        exp_data = read_file_safe(exp_file.strip("\n"))
        experimental_intensities = exp_data[:, 1]
        if magdiv:
            weights = out_file[index, : last_weight_index - 1]
            intensity = make_positive_profile(
                weights, simulated_intensities, experimental_intensities
            )
            diff_intensity = experimental_intensities - intensity
            combined_intensity.append(diff_intensity)
            # TODO: In principle these can be rewighted simulated errors
            combined_errors.append(0.01 * diff_intensity)
        else:
            diff_weights = out_file[:, last_weight_index - 1]
            print("Diff weights", diff_weights[index])
            # There maybe a case when none of this are
            if diff_weights[index] > 0.01:
                # weights = out_file[index, :last_weight_index - 1]
                # intensity = make_positive_profile(weights, simulated_intensities, experimental_intensities)
                # diff_intensity = experimental_intensities - intensity

                # if np.any(diff_simulated < 0.0):
                #    diff_simulated = shift_intensity(diff_simulated)
                intensity = diff_simulated / diff_weights[index]
                combined_intensity.append(intensity)
                combined_errors.append(0.01 * intensity)
    combined_intensity = np.average(combined_intensity, axis=0)
    combined_errors = np.average(combined_errors, axis=0)
    np.savetxt(
        os.path.join("iteration_" + str(iter), "basis_profile.dat"),
        np.transpose([qvector, combined_intensity, combined_errors]),
    )


def make_positive_profile(weights, simulated_intensities, experimental_intensity):
    """
    Correcting weights for the profile based on the inferred alpha parameter
    aplha(q) > (Isim-Iexp)/sum I_i. Final alpha parameter max[alpha(q)]
    :return:
    """
    weighted_intensities = np.dot(weights, np.transpose(simulated_intensities))
    alpha_q = (weighted_intensities - experimental_intensity) / np.sum(
        simulated_intensities
    )
    alpha_max = np.max(alpha_q)
    print("Alphas size: ", np.shape(alpha_q))
    updated_weights = weights - alpha_max
    positive_intensities = np.dot(updated_weights, np.transpose(simulated_intensities))
    return positive_intensities


def shift_intensity(intensities):
    """
    Shifting intensity so they all become positive
    :param intensities:
    :return:
    """
    min_int = np.min(intensities)
    return intensities + np.fabs(min_int)


def simple_svd(experimental):
    # Transposing due to
    intensities = experimental[:, 1].T
    errors = experimental[:, 2].T
    err_avg = np.average(errors)
    svd_a = intensities / err_avg
    svd_U, svd_s, svd_Vt = np.linalg.svd(svd_a, full_matrices=True)
    return svd_U, svd_s, svd_Vt


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
    stan_dat = {
        "sim_curves": simulated,
        "target_curve": experimental[:, 1],
        "target_errors": experimental[:, 2],
        "n_measures": np.shape(experimental)[0],
        "n_structures": np.shape(simulated)[1],
        "priors": priors,
    }

    # sm = pystan.StanModel(model_code=stan_code+psisloo_quanities)
    sm = pystan.StanModel(model_code=stan_code)
    fit = sm.sampling(
        data=stan_dat,
        iter=iterations,
        chains=chains,
        n_jobs=njobs,
        sample_file="saved_samples.txt",
    )

    # initial_values = [{"weight[0]":0.05, "weight[1]":0.1, "weight[2]":0.15,
    #                   "weight[3]":0.3, "weight[4]":0.4, "scale":1}]
    # fit = sm.optimizing(data=stan_dat, init=initial_values, algorithm="BFGS")

    # fig = fit.plot(pars="weights")
    # ax.set_color_cycle(['red', 'black', 'yellow', 'green', 'blue'])
    # fig.subplots_adjust(wspace=0.8)
    # fig.savefig("stan_weights.png", dpi=300)

    # fig = fit.plot(pars="scale")
    # fig.subplots_adjust(wspace=0.8)
    # fig.savefig("stan_scale.png", dpi=300)

    # np.savetxt("target_curve_full.csv", fit.summary()['summary'][-869:-1][:,:2])
    return fit


def calculate_stats(
    fit, experimental, simulated, cs_simulated=None, cs_rms=None, cs_experimental=None
):
    """
    Calculates statistics based on stan model
    :param fit:
    :return:
    """
    # la = fit.extract(permuted=True)  # return a dictionary of arrays
    # mu = la['weights']

    ## return an array of three dimensions: iterations, chains, parameters
    results_array = fit.extract(permuted=False, inc_warmup=False)

    nsamples = 0
    jsd_sum = 0.0
    bayesian_weights = np.zeros(np.shape(simulated)[1])
    for iteration in results_array:
        for parameters in iteration:
            current_weights = parameters[: np.shape(simulated)[1]]
            bayesian_weights += current_weights
            nsamples += 1
    bayesian_weights = bayesian_weights / nsamples

    for iteration in results_array:
        for parameters in iteration:
            current_weights = parameters[: np.shape(simulated)[1]]
            jsd_sum += JensenShannonDiv(current_weights, bayesian_weights)
    jsd = np.sqrt(jsd_sum / nsamples)

    crysol_chi2 = calculateChiCrysol(
        np.dot(bayesian_weights, np.transpose(simulated)),
        experimental[:, 1],
        experimental[:, 2],
    )
    try:
        if cs_experimental.any() != None:
            chemical_shifts_on = True
    except:
        chemical_shifts_on = False

    # scale = fit.summary(pars='scale')['summary'][0][0]
    # print("Scale from summary", scale)
    # combine_curve(experimental, simulated, bayesian_weights, scale)

    if chemical_shifts_on:
        chemshift_chi2 = calculateChemShiftsChi(
            np.dot(bayesian_weights, np.transpose(cs_simulated)),
            cs_experimental[:, 0],
            cs_experimental[:, 1],
            cs_rms,
        )
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


def execute(experimental, simulated, priors, iterations, chains, njobs):
    fit = execute_stan(experimental, simulated, priors, iterations, chains, njobs)

    return fit


def calculatePostChi():
    # Calculate chi2 between reference spectra and obtained basis spectra
    if options.diff_spectra:
        basis_data = read_file_safe(os.path.join(iteration_folder, "basis_profile.dat"))
        reference_data = read_file_safe("reference_10mer.dat")
        basis_int = basis_data[:, 1]
        reference_int = reference_data[:, 1]
        reference_err = reference_data[:, 2]
        chi2 = calculateChiCrysol(basis_int, reference_int, reference_err)
        print("chi2 with reference: ", chi2)


if __name__ == "__main__":
    doc = """
        Python interface to Complete Bayesian algorithm
        Usage: python runVBW.py --help
    """
    print(doc)
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class(usage=usage, version="0.1")

    parser.add_option(
        "-p",
        "--priors",
        dest="priors_file",
        default=None,
        help="Prior weights [OBLIGATORY]",
    )
    parser.add_option(
        "-s",
        "--simulated",
        dest="simulated_file",
        help="Simulated SAXS curves [OBLIGATORY]",
    )
    parser.add_option(
        "-S",
        "--ctd_simulated",
        dest="ctd_simulated_file",
        help="Simulated SAXS curves for CTD domain[OBLIGATORY]",
    )
    parser.add_option(
        "-e",
        "--experimental",
        dest="experimental_files",
        help="Experimental SAXS curves [OBLIGATORY]",
    )
    parser.add_option(
        "-E",
        "--simulated_errors",
        dest="simulated_errors_file",
        help="Simulated SAXS errors [OBLIGATORY]",
    )
    parser.add_option("-o", "--output", dest="output", help="Output file [OBLIGATORY]")
    parser.add_option(
        "-f",
        "--file_names",
        dest="names_file",
        default=None,
        help="File containing names of structures to process",
    )
    parser.add_option(
        "-i",
        "--iterations",
        dest="iterations",
        default=2000,
        type="int",
        help="Number of iterations",
    )
    parser.add_option(
        "-j",
        "--jobs",
        dest="njobs",
        default=8,
        type="int",
        help="Number of procceses to run",
    )
    parser.add_option(
        "-c", "--chains", dest="chains", default=4, type="int", help="Number of chains"
    )
    parser.add_option(
        "-d",
        "--diff",
        dest="diff_spectra",
        default=0,
        type="int",
        help="Wherther eneable fit of difference spectra",
    )
    parser.add_option(
        "-I",
        "--iteration",
        dest="iteration",
        default=0,
        type="int",
        help="Iteration number",
    )
    parser.add_option(
        "-a",
        "--amplitude",
        dest="amplitude",
        default=0,
        type="int",
        help="use difference spectra amplitude",
    )
    options, args = parser.parse_args()

    njobs = options.njobs
    iterations = options.iterations
    chains = options.chains
    simulated = read_file_safe(options.simulated_file)
    simulated_errors = read_file_safe(options.simulated_errors_file)
    number_of_states = np.shape(simulated)[1]
    # ctd_simulated = read_file_safe(options.ctd_simulated_file)
    file_names = read_file_safe(options.names_file, "unicode")
    # First flat prior then updated with values from previous iteration
    priors = read_file_safe(options.priors_file)
    # Initialize stan model
    if options.diff_spectra:
        # Even distribution of priors
        priors = np.ones(number_of_states + 1) * 1.0 / (number_of_states + 1)
        file_names = np.append(file_names, "diff")
    # Pickling for teh first time compilation
    sm = pystan.StanModel(model_code=stan_code)
    # with open('model.pkl', 'wb') as f:
    #    pickle.dump(sm, f)
    # sm = pickle.load(open('model.pkl', 'rb'))

    iteration_folder = "iteration_" + str(options.iteration)
    # previous_iteration_folder = 'iteration_' + str(options.iteration-1)
    if os.path.exists(iteration_folder) == False:
        os.system("mkdir " + iteration_folder)
    fout = open(os.path.join(iteration_folder, options.output), "w")
    # for name in file_names:
    #    fout.write(name + ', ')
    # fout.write('jsd, chi2\n')
    # List of experiemntal files to scan over

    stan_initialized = False

    # TODO: For each frame load diff file and combine it with simulated
    # What to do with CTD simulated?
    # make if statement for non-diff cases
    if options.diff_spectra:
        # Division by signal magnitude in the first round from second will turn to false
        average_basis_spectra(
            open(options.experimental_files).readlines(),
            simulated,
            number_of_states,
            options.iteration,
            options.amplitude,
        )
        basis_file_name = os.path.join(iteration_folder, "basis_profile.dat")
        basis_file = read_file_safe(basis_file_name)
        basis_simulated = basis_file[:, 1]
        basis_error = basis_file[:, 2]
    experimental_files = open(options.experimental_files).readlines()
    for index, experimental_file in enumerate(experimental_files):
        print("===> Analyzing frame", experimental_file)
        experimental = read_file_safe(experimental_file.strip("\n"))
        qvector = experimental[:, 0]
        errors = np.c_[simulated_errors, experimental[:, 2]]

        # threshold = 550
        # experimental = experimental[:threshold]
        # simulated = simulated[:threshold]
        # qvector = qvector[:threshold]
        # ctd_simulated = ctd_simulated[:threshold]

        target_errors = np.sqrt(np.sum(np.power(errors, 2), axis=1))
        full_simulated = simulated

        if options.diff_spectra:
            # sqrt of squared elements
            errors = np.c_[simulated_errors, basis_error, experimental[:, 2]]
            full_simulated = np.c_[simulated, basis_simulated]
            target_errors = np.sqrt(np.sum(np.power(errors, 2), axis=1))
        if priors.ndim > 1:
            print(np.shape(priors))
            priors_i = priors[index, :]
            print(priors_i)
        else:
            priors_i = priors
        # TODO: Priors should be defined as a matrix and defined for specific matrix
        stan_dat = {
            "sim_curves": full_simulated,
            # "ctd_curves" : ctd_simulated,
            "qvector": qvector,
            "target_curve": experimental[:, 1],
            "target_errors": target_errors,
            "n_measures": np.shape(experimental)[0],
            "n_structures": np.shape(full_simulated)[1],
            "priors": priors_i,
        }

        fit = sm.sampling(
            data=stan_dat, iter=iterations, chains=chains, n_jobs=njobs, refresh=-1
        )

        print(fit)
        # Calculate core shell and concatenate to simulated matrix
        # radius_shell = fit.summary(pars='radius_shell')['summary'][0][0]
        # radius_core = fit.summary(pars='radius_core')['summary'][0][0]
        # Iqcs = core_shell_model(qvector, radius_core, radius_shell)
        # print(Iqcs[0])
        # simulated_ext = np.c_[simulated, Iqcs]

        simulated_ext = full_simulated
        bayesian_weights, jsd, crysol_chi2 = calculate_stats(
            fit, experimental, simulated_ext
        )

        weights_ext = bayesian_weights
        if options.diff_spectra:
            weights_ext = bayesian_weights[:-1]
            simulated_ext = full_simulated[:, 0:-1]

        combine_curve(
            experimental,
            simulated_ext,
            weights_ext,
            options.iteration,
            experimental_file.split(".dat")[0] + ".fit",
        )

        # TODO: Write uncertainities here probably to a separate file
        for index, fname in enumerate(file_names):
            fout.write(str(bayesian_weights[index]) + " ")
        fout.write(str(jsd) + " ")
        fout.write(str(crysol_chi2) + "\n")
        # priors = bayesian_weights
        # And soem stan utilitry stats
        stan_utility.check_treedepth(fit)
        stan_utility.check_energy(fit)
        stan_utility.check_div(fit)
    fout.close()
    calculatePostChi()
    # plot_intermediates(options.output)
    # plot_chi2(options.output)
    # #Simple svd for the comparison
    # for experimental_file in open(options.experimental_files).readlines():
    #     experimental = read_file_safe(experimental_file.strip("\n"))
    #     svd_U, svd_s, svd_Vt = simple_svd(experimental)
    #     print(svd_U, svd_s, svd_Vt)
