import os
import sys
import numpy as np
import shlex
from subprocess import Popen, PIPE

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

def process_pdbs_with_experimental(pdb_list, experimental_file):
    intensities = None
    errors = None
    for pdb_file in pdb_list:
        cmd_line = 'Pepsi-SAXS --dp_min 0 --dp_max 0 --dp_N 1 --r0_min_factor 1 --r0_max_factor 1 --r0_N 1 '+\
                   pdb_file.strip('\n')+' '+experimental_file
        print(cmd_line)
        proc = Popen(shlex.split(cmd_line), stdout=PIPE, shell=False)
        os.waitpid(proc.pid, 0)
        (out, err) = proc.communicate()
        data = read_file_safe(pdb_file[:-5]+'-'+experimental_file[:-3]+'fit', 6)
        intensity = data[:,3]
        qvector = data[:,0]
        log_file = open(pdb_file[:-5]+'-'+experimental_file[:-3]+'log')
        log_data = log_file.readlines()[62]
        log_file.close()
        if log_data[:7] != 'Scaling':
            break
        scaling_factor = float(log_data.split(":")[1])
        intensity /= scaling_factor
        if intensities is None:
            intensities = intensity
        else:
            intensities = np.c_[intensities, intensity]
    return intensities

def add_errors_m(intensities, qvector):
    k = 4700
    c = 0.85
    q_arb = 0.2
    mult_fact = 100
    tau = 1
    I_0 = intensities[0]
    I_Q = intensities*mult_fact/I_0
    sigma = np.sqrt((1 / (k * qvector)) * (I_Q) * tau)
    idx = (np.abs(qvector - q_arb)).argmin()
    addon = np.sqrt((1 / (k * qvector[idx])) * 2 * c * I_Q[idx] / (1 - c)) * tau
    sigma += addon
    #I_Q_E = np.rnorm(len(qvector),I_Q*I_0/mult_fact,sigma*I_0/mult_fact)
    return sigma * I_0 / mult_fact

def generate_file_list(pdb_list):
    files = read_file_safe(pdb_list,0,'unicode')
    f = open('file_list', 'w')
    line_to_write = ' '.join(files)
    f.writelines(line_to_write)
    f.close()

def generate_weights(pdb_list):
    out = open("weights.txt", "w")
    files = read_file_safe(pdb_list, 0, 'unicode')
    len_models = len(files)
    fwght = 1.0 / (len_models)
    for i in range(len_models):
        out.write(str(fwght) + " ")
    out.close()

if __name__ == "__main__":
    pdb_list_name = sys.argv[1]
    pdb_list = open(pdb_list_name).readlines()
    experimental_file = sys.argv[2]
    generate_file_list(pdb_list_name)
    generate_weights(pdb_list_name)
    intensities = process_pdbs_with_experimental(pdb_list, experimental_file)
    np.savetxt('SimulatedIntensities.txt', intensities)