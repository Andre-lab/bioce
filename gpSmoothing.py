"""
Smoothing curve using Gaussian process.
Taking output from VBI as a input and process it.
Returns sets of new weights for each time point.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, ExpSineSquared, RationalQuadratic, ConstantKernel, Matern


def read_file_safe(filename, skip_lines, dtype="float64"):
    """
    Simple check if file exists
    :param filename:
    :return:
    """
    try:
        results = np.genfromtxt(filename, skip_header=skip_lines, dtype=dtype)
    except IOError as err:
        print(os.strerror(err.errno))
    return results

def smooth_curve(weights, kernel):
    """
    Curve smoothing using Gaussian Process
    :param data:
    :return:
    """

    # tvector = np.array([0.05, 0.08, 0.12, 0.15, 0.22, 0.30, 0.43, 0.6, 0.92, 1.3, 1.9, 2.7, 3.96,
    #                     5.81, 8.15, 12.47, 18.30, 26.90, 39.45, 57.95])
    tvector = np.linspace(1, 20, 20)
    xshape = np.shape(weights)[0]
    x = np.atleast_2d(np.linspace(0, 21, 1000)).T
    X = np.reshape(tvector, (xshape, 1))
    #y = np.reshape(weights, (xshape, 1))
    y = weights
    print(np.shape(X), np.shape(y), np.shape(x))
    gp = GaussianProcessRegressor(kernel=kernel)
    # smoothed_curve, sigma = gp.predict(X, return_std=True)
    # Fit to data using Maximum Likelihood Estimation of the parameters
    gp.fit(X, y)
    # Make the prediction on the meshed x-axis (ask for MSE as well)
    smoothed_curve, sigma = gp.predict(x, return_std=True)
    print(np.shape(smoothed_curve), np.shape(sigma))
    y_samples = gp.sample_y(x, 10)
    return tvector, smoothed_curve, sigma, y_samples

def plot_curve(tvector, weights, smoothed_weigths, sigma, frame_index, y_samples):
    """
    Produces plot given the data
    :param data:
    :param log:
    :return:
    """
    # print(np.shape(tvector))
    # print(np.shape(weights))
    # print(np.shape(smoothed_weights))
    # print(np.shape(sigma))

    matplotlib.use("TkAgg")
    plt.figure()
    x = np.atleast_2d(np.linspace(0, 21, 1000)).T
    plt.plot(tvector, weights, 'ro', markersize=4, mfc="none", label="Inferred weights")
    plt.plot(x, smoothed_weigths, 'b-', label="Smoothed weights")
    #plt.fill_between(tvector, smoothed_weigths - sigma, smoothed_weigths + sigma,
    #                 alpha=0.2, color='k')
    plt.fill(np.concatenate([x, x[::-1]]),
                   np.concatenate([smoothed_weigths - 1.9600 * sigma,
                                   (smoothed_weigths + 1.9600 * sigma)[::-1]]),
                   alpha=.5, fc='b', ec='None', label='95% confidence interval')


    #plt.plot(tvector, y_samples, lw=1)

    plt.ylabel("$weights$")
    plt.xlabel("$timeframe$")
    plt.legend(loc='upper right')
    plt.ylim(-0.1, 0.5)
    plt.savefig("smooth_curve_"+str(index)+".png", dpi=300, bbox_inches='tight')
    if frame_index == 0:
        plt.show()


if __name__ == "__main__":

    kernels = [1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-1, 10.0)),
               1.0 * RationalQuadratic(length_scale=1.0, alpha=0.1),
               1.0 * ExpSineSquared(length_scale=1.0, periodicity=3.0,
                                    length_scale_bounds=(0.1, 10.0),
                                    periodicity_bounds=(1.0, 10.0)),
               ConstantKernel(0.1, (0.01, 10.0))
               * (DotProduct(sigma_0=1.0, sigma_0_bounds=(0.1, 10.0)) ** 2),
               1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0),
                            nu=1.5),
               DotProduct() + WhiteKernel(1e-4)]

    data = read_file_safe('batch_test',0)
    #Read in and remove last column for chi2, then run smoothing in the loop
    number_of_species = np.shape(data)[1] - 2
    #number_of_species = 1
    for index in range(number_of_species):
        weights = data[:,index]
        #for kernel in kernels:
        #kernel = ConstantKernel(0.002, (0.002, 1.0))* (DotProduct(sigma_0=1.0, sigma_0_bounds=(0.1, 1.0)) ** 2)
        #kernel = DotProduct() + WhiteKernel(1e-5)
        kernel = 1.0 * RBF(length_scale=0.7, length_scale_bounds=(1e-5, 1e5)) \
              + WhiteKernel(noise_level=1e-3, noise_level_bounds=(1e-5, 1e+5))
        #kernel = 1.0 * ExpSineSquared(length_scale=1, periodicity=4.0,
        #                      length_scale_bounds=(0.1, 10.0),
        #                      periodicity_bounds=(1.0, 10.0)) * ConstantKernel(0.1, (0.01, 10.0))
        try:
            tvector, smoothed_weights, sigma, y_samples = smooth_curve(weights, kernel)
            plot_curve(tvector, weights, smoothed_weights, sigma, index, y_samples)
            smoothed_weights[smoothed_weights < 0.0] = 0.01
            if index == 0:
                new_priors = smoothed_weights
            else:
                new_priors = np.vstack([new_priors, smoothed_weights])
        except:
            pass
    new_priors = np.transpose(new_priors)
    row_sums = new_priors.sum(axis=1)
    new_priors = new_priors / row_sums[:, np.newaxis]
    # Scaling factors for alphas
    scaling_factor = 10
    new_priors = scaling_factor * new_priors
    np.savetxt('prior_matrix.txt',new_priors)

    #test_sk()