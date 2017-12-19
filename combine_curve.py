import numpy as np
import optparse
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

def combine_curve(simulated, weights, scale):
    """

    :param simulated:
    :param weights:
    :param scale:
    :return:
    """
    combuined = scale*np.dot(weights, np.transpose(simulated))

    return combuined

if __name__ == "__main__":
    doc = """
        Python interface to Variational Bayesian algorithm
        Usage: python runVBW.py --help
    """
    print(doc)
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-w", "--weights", dest="weights", default=None,
                      help="Posterior weights [OBLIGATORY]")
    parser.add_option("-s", "--simulated", dest="simulated_file",
                      help="Simulated SAXS curves [OBLIGATORY]")
    parser.add_option("-e", "--experimenttal", dest="experimental_file",
                      help="Simulated SAXS curves [OBLIGATORY]")
    parser.add_option("-k", "--scale", dest="scale", default=1,
                  type='float',
                  help="scaling factor")
    options, args = parser.parse_args()
    scale = options.scale
    weights = read_file_safe(options.weights)
    simulated = read_file_safe(options.simulated_file)
    experimental = read_file_safe(options.experimental_file)
    combined = combine_curve(simulated, weights, scale)
    q_column = experimental[:,0]
    np.savetxt("combinedCurveqIq.txt", np.transpose((q_column, combined)))